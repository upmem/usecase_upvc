/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <numa.h>

#include "accumulateread.h"
#include "dispatch.h"
#include "dpu_backend.h"
#include "genome.h"
#include "getread.h"
#include "index.h"
#include "parse_args.h"
#include "processread.h"
#include "simu_backend.h"
#include "upvc.h"
#include "vartree.h"

#include "backends_functions.h"

FILE *time_file;
pthread_mutex_t time_file_mutex;

unsigned int nb_ranks_per_run;
unsigned int nb_dpus_per_run;

static backends_functions_t backends_functions;
static unsigned int round;
static FILE *fipe1, *fipe2, *fope1, *fope2;
static sem_t *getreads_to_dispatch_sem, *dispatch_to_exec_sem, *exec_to_dispatch_sem, *exec_to_acc_sem,
    *accprocess_to_getreads_sem, *acc_to_process_sem;
sem_t *acc_to_exec_sem;

#define DEBUG_DPU_OFFSET 2412
#define DEBUG_PASS 3

#define LAST_RUN(dpu_offset) (((dpu_offset) + nb_dpus_per_run) >= index_get_nb_dpu())
#ifndef DEBUG_DPU_OFFSET
#define FOREACH_RUN(dpu_offset) for (unsigned int dpu_offset = 0; dpu_offset < index_get_nb_dpu(); dpu_offset += nb_dpus_per_run)
#else
#define FOREACH_RUN(dpu_offset) __attribute__((unused)) unsigned int dpu_offset = DEBUG_DPU_OFFSET;
#endif
#define FOREACH_RANK(each_rank) for (unsigned int each_rank = 0; each_rank < nb_ranks_per_run; each_rank++)
#ifndef DEBUG_PASS
#define FOREACH_PASS(each_pass) for (unsigned int each_pass = 0; get_reads_in_buffer(each_pass) != 0; each_pass++)
#else
#define FOREACH_PASS(each_pass) unsigned int each_pass = DEBUG_PASS;
#endif
#define FOR(loop) for (int _i = 0; _i < (loop); _i++)

void *thread_get_reads(__attribute__((unused)) void *arg)
{
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_READS_BUFFER) { sem_post(accprocess_to_getreads_sem); }

        unsigned int each_pass = 0;
        do {
#ifdef DEBUG_PASS
            for (each_pass = 0; each_pass < DEBUG_PASS; each_pass++) {
                get_reads(fipe1, fipe2, each_pass);
            }
#endif
            sem_wait(accprocess_to_getreads_sem);

            PRINT_TIME_GET_READS();
            get_reads(fipe1, fipe2, each_pass);
            PRINT_TIME_GET_READS();

            sem_post(getreads_to_dispatch_sem);
#ifdef DEBUG_PASS
            return NULL;
#endif
        } while (get_reads_in_buffer(each_pass++) != 0);

        fseek(fipe1, 0, SEEK_SET);
        fseek(fipe2, 0, SEEK_SET);

        FOR(NB_READS_BUFFER) { sem_wait(accprocess_to_getreads_sem); }
    }

    return NULL;
}

void *thread_exec_rank(void *arg)
{
    double t1, t2;
    const unsigned int rank_id = *(unsigned int *)arg;
    const unsigned int delta_neighbour = (SIZE_SEED * round) / 4;

    //int ret = numa_run_on_node(3 - get_rank_numa_node(rank_id));
    int ret = numa_run_on_node(get_rank_numa_node(rank_id));
    if (ret < 0) {
        printf("ERROR numa failed\n");
        return NULL;
    }

    FOREACH_RUN(dpu_offset)
    {
        print_line(rank_id);
        print(rank_id, "M %u", dpu_offset);
        PRINT_TIME_WRITE_MRAM(rank_id);
        t1 = my_clock();
        backends_functions.load_mram(dpu_offset, rank_id, delta_neighbour);
        t2 = my_clock();
        PRINT_TIME_WRITE_MRAM(rank_id);
        print(rank_id, "%.2lf", t2 - t1);

        sem_wait(&dispatch_to_exec_sem[rank_id]);

        FOREACH_PASS(each_pass)
        {
            print_line(rank_id);

            print(rank_id, "P %u", each_pass);
            PRINT_TIME_MAP_READ(rank_id);
            t1 = my_clock();
            backends_functions.run_dpu(dpu_offset, rank_id, each_pass, &exec_to_dispatch_sem[rank_id], &acc_to_exec_sem[rank_id]);
            t2 = my_clock();
            PRINT_TIME_MAP_READ(rank_id);
            print(rank_id, "T %.2lf", t2 - t1);

            sem_post(&exec_to_acc_sem[rank_id]);
            sem_wait(&dispatch_to_exec_sem[rank_id]);
        }

        sem_post(&exec_to_acc_sem[rank_id]);
    }
    return NULL;
}

void *thread_dispatch(__attribute__((unused)) void *arg)
{
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            FOREACH_RANK(each_rank) { sem_post(&exec_to_dispatch_sem[each_rank]); }
        }

        sem_wait(getreads_to_dispatch_sem);
        FOREACH_PASS(each_pass)
        {
            FOREACH_RANK(each_rank) { sem_wait(&exec_to_dispatch_sem[each_rank]); }

            PRINT_TIME_DISPATCH();
            dispatch_read(each_pass);
            PRINT_TIME_DISPATCH();

            FOREACH_RANK(each_rank) { sem_post(&dispatch_to_exec_sem[each_rank]); }

            sem_wait(getreads_to_dispatch_sem);
        }

        FOREACH_RANK(each_rank) { sem_post(&dispatch_to_exec_sem[each_rank]); }

        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            FOREACH_RANK(each_rank) { sem_wait(&exec_to_dispatch_sem[each_rank]); }
        }
    }
    return NULL;
}

void *thread_acc(__attribute__((unused)) void *arg)
{
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            FOREACH_RANK(each_rank) { sem_post(&acc_to_exec_sem[each_rank]); }
        }

        FOREACH_RANK(each_rank) { sem_wait(&exec_to_acc_sem[each_rank]); }

        FOREACH_PASS(each_pass)
        {
            PRINT_TIME_ACC_READ();
            accumulate_read(each_pass, dpu_offset);
            PRINT_TIME_ACC_READ();

            FOREACH_RANK(each_rank) { sem_post(&acc_to_exec_sem[each_rank]); }
            if (LAST_RUN(dpu_offset)) {
                sem_post(acc_to_process_sem);
            } else {
                sem_post(accprocess_to_getreads_sem);
            }
            FOREACH_RANK(each_rank) { sem_wait(&exec_to_acc_sem[each_rank]); }
        }
        if (LAST_RUN(dpu_offset)) {
            sem_post(acc_to_process_sem);
        } else {
            sem_post(accprocess_to_getreads_sem);
        }
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            FOREACH_RANK(each_rank) { sem_wait(&acc_to_exec_sem[each_rank]); }
        }
    }
    return NULL;
}

void *thread_process(__attribute__((unused)) void *arg)
{
    sem_wait(acc_to_process_sem);

    FOREACH_PASS(each_pass)
    {
        PRINT_TIME_PROCESS_READ();
        process_read(fope1, fope2, round, each_pass);
        PRINT_TIME_PROCESS_READ();
        sem_post(accprocess_to_getreads_sem);

        sem_wait(acc_to_process_sem);
    }

    sem_post(accprocess_to_getreads_sem);

    return NULL;
}

static void exec_round()
{
    char filename[FILENAME_MAX];
    char *input_prefix = get_input_path();
    static unsigned int max_nb_pass;

    if (round == 0) {
        size_t read_size1, read_size2, nb_read1, nb_read2;

        sprintf(filename, "%s_PE1.fastq", input_prefix);
        fipe1 = fopen(filename, "r");
        CHECK_FILE(fipe1, filename);
        assert(get_input_info(fipe1, &read_size1, &nb_read1) == 0);
        assert(read_size1 == SIZE_READ);

        sprintf(filename, "%s_PE2.fastq", input_prefix);
        fipe2 = fopen(filename, "r");
        CHECK_FILE(fipe2, filename);
        assert(get_input_info(fipe2, &read_size2, &nb_read2) == 0);
        assert(read_size2 == SIZE_READ);

        assert(nb_read1 == nb_read2);
        max_nb_pass = (unsigned int)(nb_read1 * 2 + nb_read2 * 2 + MAX_READS_BUFFER - 1) / MAX_READS_BUFFER;
    } else {
        fipe1 = fope1;
        fipe2 = fope2;
    }

    if (round == NB_ROUND - 1) {
        fope1 = NULL;
        fope2 = NULL;
    } else {
        sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round + 1);
        fope1 = fopen(filename, "w+");
        assert(fope1 != NULL);
        assert(unlink(filename) == 0);
        sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round + 1);
        fope2 = fopen(filename, "w+");
        assert(fope2 != NULL);
        assert(unlink(filename) == 0);
    }

    accumulate_init(max_nb_pass);

    pthread_t tid_get_reads;
    pthread_t tid_exec_rank[nb_ranks_per_run];
    pthread_t tid_dispatch;
    pthread_t tid_acc;
    pthread_t tid_process;

    unsigned int exec_rank_arg[nb_ranks_per_run];
    FOREACH_RANK(each_rank) { exec_rank_arg[each_rank] = each_rank; }

    int ret;

    // INIT
    getreads_to_dispatch_sem = malloc(sizeof(sem_t));
    dispatch_to_exec_sem = malloc(sizeof(sem_t) * nb_ranks_per_run);
    exec_to_dispatch_sem = malloc(sizeof(sem_t) * nb_ranks_per_run);
    exec_to_acc_sem = malloc(sizeof(sem_t) * nb_ranks_per_run);
    acc_to_exec_sem = malloc(sizeof(sem_t) * nb_ranks_per_run);
    accprocess_to_getreads_sem = malloc(sizeof(sem_t));
    acc_to_process_sem = malloc(sizeof(sem_t));
    assert(getreads_to_dispatch_sem != NULL && dispatch_to_exec_sem != NULL && exec_to_dispatch_sem != NULL
        && exec_to_acc_sem != NULL && acc_to_exec_sem != NULL && accprocess_to_getreads_sem != NULL
        && acc_to_process_sem != NULL);

    ret = sem_init(getreads_to_dispatch_sem, 0, 0);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = sem_init(&dispatch_to_exec_sem[each_rank], 0, 0);
        assert(ret == 0);
        ret = sem_init(&exec_to_dispatch_sem[each_rank], 0, 0);
        assert(ret == 0);
        ret = sem_init(&exec_to_acc_sem[each_rank], 0, 0);
        assert(ret == 0);
        ret = sem_init(&acc_to_exec_sem[each_rank], 0, 0);
        assert(ret == 0);
    }
    ret = sem_init(acc_to_process_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(accprocess_to_getreads_sem, 0, 0);
    assert(ret == 0);

    // CREATE
    ret = pthread_create(&tid_get_reads, NULL, thread_get_reads, NULL);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = pthread_create(&tid_exec_rank[each_rank], NULL, thread_exec_rank, (void *)&exec_rank_arg[each_rank]);
        assert(ret == 0);
    }
    ret = pthread_create(&tid_dispatch, NULL, thread_dispatch, NULL);
    assert(ret == 0);
    ret = pthread_create(&tid_acc, NULL, thread_acc, NULL);
    assert(ret == 0);
    ret = pthread_create(&tid_process, NULL, thread_process, NULL);
    assert(ret == 0);

    // JOIN
    ret = pthread_join(tid_get_reads, NULL);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = pthread_join(tid_exec_rank[each_rank], NULL);
        assert(ret == 0);
    }
    ret = pthread_join(tid_dispatch, NULL);
    assert(ret == 0);
    ret = pthread_join(tid_acc, NULL);
    assert(ret == 0);
    ret = pthread_join(tid_process, NULL);
    assert(ret == 0);

    // DESTROY
    ret = sem_destroy(getreads_to_dispatch_sem);
    free(getreads_to_dispatch_sem);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = sem_destroy(&dispatch_to_exec_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&exec_to_dispatch_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&acc_to_exec_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&exec_to_acc_sem[each_rank]);
        assert(ret == 0);
    }
    free(dispatch_to_exec_sem);
    free(exec_to_dispatch_sem);
    free(acc_to_exec_sem);
    free(exec_to_acc_sem);
    ret = sem_destroy(acc_to_process_sem);
    free(acc_to_process_sem);
    assert(ret == 0);
    ret = sem_destroy(accprocess_to_getreads_sem);
    free(accprocess_to_getreads_sem);
    assert(ret == 0);

    accumulate_free();

    fclose(fipe1);
    fclose(fipe2);
}

static void init_time_file_and_mutex()
{
    char filename[1024];
    sprintf(filename, "%s_time.csv", get_input_path());
    time_file = fopen(filename, "w");
    CHECK_FILE(time_file, filename);
    fprintf(time_file,
        "time, get_reads, dispatch, accumulate_read, process_read, "
        "write_mram, write_reads, compute, read_result, map_read\n");

    pthread_mutex_init(&time_file_mutex, NULL);
}

static void close_time_file_and_mutex()
{
    fclose(time_file);
    pthread_mutex_destroy(&time_file_mutex);
}

static void do_mapping()
{
    variant_tree_init();
    init_time_file_and_mutex();
    backends_functions.init_backend(&nb_dpus_per_run, &nb_ranks_per_run);
    dispatch_init();
    process_read_init();

    for (round = 0; round < NB_ROUND; round++) {
        printf("#################\n"
               "starting round %u\n"
               "#################\n",
            round);
        exec_round();
    }
    create_vcf();

    process_read_free();
    dispatch_free();
    backends_functions.free_backend();
    close_time_file_and_mutex();
    variant_tree_free();
}

static void print_time()
{
    time_t timer;
    char time_buf[26];
    struct tm *tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    printf("upvc started at: %s\n", time_buf);
}

int main(int argc, char *argv[])
{
    validate_args(argc, argv);

    printf("%s\n", VERSION);
    print_time();

    printf("Information:\n");
    printf("\tread size: %d\n", SIZE_READ);
    printf("\tseed size: %u\n", SIZE_SEED);

    if (get_simulation_mode()) {
        backends_functions.init_backend = init_backend_simulation;
        backends_functions.free_backend = free_backend_simulation;
        backends_functions.run_dpu = run_dpu_simulation;
        backends_functions.load_mram = load_mram_simulation;
    } else {
        backends_functions.init_backend = init_backend_dpu;
        backends_functions.free_backend = free_backend_dpu;
        backends_functions.run_dpu = run_on_dpu;
        backends_functions.load_mram = load_mram_dpu;
    }

    switch (get_goal()) {
    case goal_index:
        index_create_folder();
        genome_create();
        index_create();
        break;
    case goal_map:
        genome_load();
        index_load();
        do_mapping();
        break;
    case goal_unknown:
    default:
        ERROR_EXIT(ERR_NO_GOAL_DEFINED, "goal has not been specified!");
    }

    index_free();
    genome_free();
    free_args();

    return 0;
}
