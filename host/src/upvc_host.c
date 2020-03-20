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

#include "accumulateread.h"
#include "dispatch.h"
#include "dpu_backend.h"
#include "dpus_mgmt.h"
#include "genome.h"
#include "getread.h"
#include "index.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "processread.h"
#include "simu_backend.h"
#include "upvc.h"
#include "upvc_dpu.h"
#include "vartree.h"
#include "vcf.h"

#include "backends_functions.h"

FILE *time_file;
pthread_mutex_t time_file_mutex;
backends_functions_t backends_functions;

#define FOREACH_RUN(dpu_offset) for (unsigned int dpu_offset = 0; dpu_offset < get_nb_dpu(); dpu_offset += get_nb_dpus_per_run())
#define FOREACH_RANK(each_rank) for (unsigned int each_rank = 0; each_rank < dpu_get_nb_ranks_per_run(); each_rank++)
#define FOREACH_PASS(each_pass) for (unsigned int each_pass = 0; get_reads_in_buffer(each_pass) != 0; each_pass++)

typedef struct {
    sem_t *dispatch_free_sem;
    sem_t *acc_process_wait_sem;
    FILE *fipe1;
    FILE *fipe2;
} thread_get_reads_arg_t;

void *thread_get_reads(void *arg)
{
    thread_get_reads_arg_t *args = (thread_get_reads_arg_t *)arg;
    FILE *fipe1 = args->fipe1;
    FILE *fipe2 = args->fipe2;
    sem_t *dispatch_free_sem = args->dispatch_free_sem;
    sem_t *acc_process_wait_sem = args->acc_process_wait_sem;

    FOREACH_RUN(dpu_offset)
    {
        unsigned int each_pass = 0;
        do {
            sem_wait(acc_process_wait_sem);

            PRINT_TIME_GET_READS();
            get_reads(fipe1, fipe2, each_pass);
            PRINT_TIME_GET_READS();

            sem_post(dispatch_free_sem);
        } while (get_reads_in_buffer(each_pass++) != 0);

        fseek(fipe1, 0, SEEK_SET);
        fseek(fipe2, 0, SEEK_SET);

        // wait for run to be complete before starting the next one
        for (unsigned int each_reads_buffer = 0; each_reads_buffer < NB_READS_BUFFER - 1; each_reads_buffer++)
            sem_wait(acc_process_wait_sem);
    }

    return NULL;
}

typedef struct {
    sem_t *dispatch_wait_sem;
    sem_t *dispatch_free_sem;
    sem_t *acc_wait_sem;
    sem_t *acc_free_sem;
    unsigned int rank_id;
    int round;
} thread_exec_rank_arg_t;

void *thread_exec_rank(void *arg)
{
    double t1, t2;
    thread_exec_rank_arg_t *args = (thread_exec_rank_arg_t *)arg;
    sem_t *dispatch_wait_sem = args->dispatch_wait_sem;
    sem_t *dispatch_free_sem = args->dispatch_free_sem;
    sem_t *acc_wait_sem = args->acc_wait_sem;
    sem_t *acc_free_sem = args->acc_free_sem;
    unsigned int rank_id = args->rank_id;
    int round = args->round;

    unsigned int delta_neighbour = (SIZE_SEED * round) / 4;

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

        // wait for get_reads_in_buffer to mean something
        sem_wait(dispatch_wait_sem);

        FOREACH_PASS(each_pass)
        {
            print_line(rank_id);

            print(rank_id, "P %u", each_pass);
            PRINT_TIME_MAP_READ(rank_id);
            t1 = my_clock();
            backends_functions.run_dpu(dpu_offset, rank_id, delta_neighbour, dispatch_free_sem, acc_wait_sem);
            t2 = my_clock();
            PRINT_TIME_MAP_READ(rank_id);
            print(rank_id, "T %.2lf", t2 - t1);

            sem_post(acc_free_sem);
            sem_wait(dispatch_wait_sem);
        }

        // let acc sees that this run is complete
        sem_post(acc_free_sem);
    }
    return NULL;
}

typedef struct {
    sem_t *get_reads_wait_sem;
    sem_t *exec_rank_wait_sem;
    sem_t *exec_rank_free_sem;
} thread_dispatch_arg_t;

void *thread_dispatch(void *arg)
{
    thread_dispatch_arg_t *args = (thread_dispatch_arg_t *)arg;
    sem_t *get_reads_wait_sem = args->get_reads_wait_sem;
    sem_t *exec_rank_wait_sem = args->exec_rank_wait_sem;
    sem_t *exec_rank_free_sem = args->exec_rank_free_sem;

    FOREACH_RUN(dpu_offset)
    {
        // wait for get_reads_in_buffer to mean something
        sem_wait(get_reads_wait_sem);

        FOREACH_PASS(each_pass)
        {
            FOREACH_RANK(each_rank) { sem_wait(&exec_rank_wait_sem[each_rank]); }

            PRINT_TIME_DISPATCH();
            dispatch_read(each_pass);
            PRINT_TIME_DISPATCH();

            FOREACH_RANK(each_rank) { sem_post(&exec_rank_free_sem[each_rank]); }

            sem_wait(get_reads_wait_sem);
        }

        // let exec_rank sees that this run is complete
        FOREACH_RANK(each_rank) { sem_post(&exec_rank_free_sem[each_rank]); }
    }
    return NULL;
}

typedef struct {
    sem_t *exec_rank_wait_sem;
    sem_t *exec_rank_free_sem;
    sem_t *process_free_sem;
    sem_t *get_reads_free_sem;
    unsigned int *result_tab_nb_read;
    dpu_result_out_t **result_tab;
} thread_acc_arg_t;

void *thread_acc(void *arg)
{
    thread_acc_arg_t *args = (thread_acc_arg_t *)arg;
    sem_t *exec_rank_wait_sem = args->exec_rank_wait_sem;
    sem_t *exec_rank_free_sem = args->exec_rank_free_sem;
    sem_t *process_free_sem = args->process_free_sem;
    sem_t *get_reads_free_sem = args->get_reads_free_sem;

    unsigned int nb_dpu = get_nb_dpu();
    unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
    unsigned int dpu_offset = 0;

    for (; dpu_offset < (nb_dpu - nb_dpus_per_run); dpu_offset += nb_dpus_per_run) {
        // allow get_reads to start the run
        for (unsigned int each_reads_buffer = 0; each_reads_buffer < NB_READS_BUFFER; each_reads_buffer++)
            sem_post(get_reads_free_sem);
        // wait for get_read_in_buffer to mean something
        FOREACH_RANK(each_rank) { sem_wait(&exec_rank_wait_sem[each_rank]); }

        FOREACH_PASS(each_pass)
        {
            PRINT_TIME_ACC_READ();
            accumulate_read(each_pass, dpu_offset);
            PRINT_TIME_ACC_READ();

            sem_post(get_reads_free_sem);
            FOREACH_RANK(each_rank) { sem_post(&exec_rank_free_sem[each_rank]); }
            FOREACH_RANK(each_rank) { sem_wait(&exec_rank_wait_sem[each_rank]); }
        }
    }
    {
        // allow get_reads to start the run
        for (unsigned int each_reads_buffer = 0; each_reads_buffer < NB_READS_BUFFER; each_reads_buffer++)
            sem_post(get_reads_free_sem);
        // wait for get_read_in_buffer to mean something
        FOREACH_RANK(each_rank) { sem_wait(&exec_rank_wait_sem[each_rank]); }

        FOREACH_PASS(each_pass)
        {
            PRINT_TIME_ACC_READ();
            accumulate_read(each_pass, dpu_offset);
            PRINT_TIME_ACC_READ();

            sem_post(process_free_sem);
            FOREACH_RANK(each_rank) { sem_post(&exec_rank_free_sem[each_rank]); }
            FOREACH_RANK(each_rank) { sem_wait(&exec_rank_wait_sem[each_rank]); }
        }
        // let process sees that this run is complete
        sem_post(process_free_sem);
    }
    return NULL;
}

typedef struct {
    sem_t *acc_wait_sem;
    sem_t *get_reads_free_sem;
    FILE *fope1;
    FILE *fope2;
    int round;
} thread_process_arg_t;

void *thread_process(void *arg)
{
    thread_process_arg_t *args = (thread_process_arg_t *)arg;
    sem_t *acc_wait_sem = args->acc_wait_sem;
    sem_t *get_reads_free_sem = args->get_reads_free_sem;
    FILE *fope1 = args->fope1;
    FILE *fope2 = args->fope2;
    int round = args->round;

    // wait for get_read_in_buffer to mean something
    sem_wait(acc_wait_sem);

    FOREACH_PASS(each_pass)
    {
        PRINT_TIME_PROCESS_READ();
        process_read(fope1, fope2, round, each_pass);
        PRINT_TIME_PROCESS_READ();
        sem_post(get_reads_free_sem);

        sem_wait(acc_wait_sem);
    }

    return NULL;
}

static void exec_round(int round)
{
    char filename[1024];
    FILE *fipe1, *fipe2, *fope1, *fope2;
    char *input_prefix = get_input_path();
    unsigned int nb_rank = dpu_get_nb_ranks_per_run();

    sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round + 1);
    fope1 = fopen(filename, "w");
    sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round + 1);
    fope2 = fopen(filename, "w");

    if (round == 0) {
        sprintf(filename, "%s_PE1.fastq", input_prefix);
        fipe1 = fopen(filename, "r");
        sprintf(filename, "%s_PE2.fastq", input_prefix);
        fipe2 = fopen(filename, "r");
    } else {
        sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round);
        fipe1 = fopen(filename, "r");
        sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round);
        fipe2 = fopen(filename, "r");
    }

    pthread_t tid_get_reads;
    pthread_t tid_exec_rank[nb_rank];
    pthread_t tid_dispatch;
    pthread_t tid_acc;
    pthread_t tid_process;

    sem_t get_reads_dispatch_sem;
    sem_t exec_rank_dispatch_sem[nb_rank];
    sem_t exec_rank_acc_sem[nb_rank];
    sem_t dispatch_exec_rank_sem[nb_rank];
    sem_t acc_exec_rank_sem[nb_rank];
    sem_t acc_process_sem;
    sem_t acc_process_get_reads_sem;

    thread_get_reads_arg_t get_reads_arg = {
        .dispatch_free_sem = &get_reads_dispatch_sem,
        .acc_process_wait_sem = &acc_process_get_reads_sem,
        .fipe1 = fipe1,
        .fipe2 = fipe2,
    };
    thread_exec_rank_arg_t exec_rank_arg[nb_rank];
    FOREACH_RANK(each_rank)
    {
        exec_rank_arg[each_rank].dispatch_wait_sem = &dispatch_exec_rank_sem[each_rank];
        exec_rank_arg[each_rank].dispatch_free_sem = &exec_rank_dispatch_sem[each_rank];
        exec_rank_arg[each_rank].acc_wait_sem = &acc_exec_rank_sem[each_rank];
        exec_rank_arg[each_rank].acc_free_sem = &exec_rank_acc_sem[each_rank];
        exec_rank_arg[each_rank].rank_id = each_rank;
        exec_rank_arg[each_rank].round = round;
    }
    thread_dispatch_arg_t dispatch_arg = {
        .get_reads_wait_sem = &get_reads_dispatch_sem,
        .exec_rank_wait_sem = exec_rank_dispatch_sem,
        .exec_rank_free_sem = dispatch_exec_rank_sem,
    };
    thread_acc_arg_t acc_arg = {
        .exec_rank_wait_sem = exec_rank_acc_sem,
        .exec_rank_free_sem = acc_exec_rank_sem,
        .process_free_sem = &acc_process_sem,
        .get_reads_free_sem = &acc_process_get_reads_sem,
    };
    thread_process_arg_t process_arg = {
        .acc_wait_sem = &acc_process_sem,
        .get_reads_free_sem = &acc_process_get_reads_sem,
        .fope1 = fope1,
        .fope2 = fope2,
        .round = round,
    };
    int ret;

    // INIT
    ret = sem_init(&get_reads_dispatch_sem, 0, 0);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = sem_init(&exec_rank_dispatch_sem[each_rank], 0, 1);
        assert(ret == 0);
        ret = sem_init(&exec_rank_acc_sem[each_rank], 0, 0);
        assert(ret == 0);
        ret = sem_init(&dispatch_exec_rank_sem[each_rank], 0, 0);
        assert(ret == 0);
        ret = sem_init(&acc_exec_rank_sem[each_rank], 0, 1);
        assert(ret == 0);
    }
    ret = sem_init(&acc_process_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&acc_process_get_reads_sem, 0, 0);
    assert(ret == 0);

    // CREATE
    ret = pthread_create(&tid_get_reads, NULL, thread_get_reads, (void *)&get_reads_arg);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = pthread_create(&tid_exec_rank[each_rank], NULL, thread_exec_rank, (void *)&exec_rank_arg[each_rank]);
        assert(ret == 0);
    }
    ret = pthread_create(&tid_dispatch, NULL, thread_dispatch, (void *)&dispatch_arg);
    assert(ret == 0);
    ret = pthread_create(&tid_acc, NULL, thread_acc, (void *)&acc_arg);
    assert(ret == 0);
    ret = pthread_create(&tid_process, NULL, thread_process, (void *)&process_arg);
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
    ret = sem_destroy(&get_reads_dispatch_sem);
    assert(ret == 0);
    FOREACH_RANK(each_rank)
    {
        ret = sem_destroy(&exec_rank_dispatch_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&exec_rank_acc_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&dispatch_exec_rank_sem[each_rank]);
        assert(ret == 0);
        ret = sem_destroy(&acc_exec_rank_sem[each_rank]);
        assert(ret == 0);
    }
    ret = sem_destroy(&acc_process_sem);
    assert(ret == 0);
    ret = sem_destroy(&acc_process_get_reads_sem);
    assert(ret == 0);

    fclose(fipe1);
    fclose(fipe2);
    fclose(fope1);
    fclose(fope2);
}

static void init_time_file_and_mutex()
{
    char filename[1024];
    sprintf(filename, "%s_time.csv", get_input_path());
    time_file = fopen(filename, "w");
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
    init_time_file_and_mutex();
    backends_functions.init_backend();
    unsigned int nb_dpu = get_nb_dpu();
    dispatch_init(nb_dpu);

    for (int round = 0; round < 3; round++) {
        printf("#################\n"
               "starting round %u\n"
               "#################\n",
            round);
        accumulate_init();
        exec_round(round);
        accumulate_free();
    }
    create_vcf();

    dispatch_free(nb_dpu);
    backends_functions.free_backend(nb_dpu);
    close_time_file_and_mutex();
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

    assert(get_read_size(get_input_pe1()) == SIZE_READ);
    printf("Information:\n");
    printf("\tread size: %d\n", SIZE_READ);

    if (get_simulation_mode()) {
        backends_functions.init_backend = init_backend_simulation;
        backends_functions.free_backend = free_backend_simulation;
        backends_functions.run_dpu = run_dpu_simulation;
        backends_functions.add_seed_to_requests = add_seed_to_simulation_requests;
        backends_functions.load_mram = load_mram_simulation;
    } else {
        backends_functions.init_backend = init_backend_dpu;
        backends_functions.free_backend = free_backend_dpu;
        backends_functions.run_dpu = run_on_dpu;
        backends_functions.add_seed_to_requests = add_seed_to_dpu_requests;
        backends_functions.load_mram = load_mram_dpu;
    }

    genome_init(get_input_fasta());

    switch (get_goal()) {
    case goal_index:
        index_init(get_nb_dpu());
        index_save();
        break;
    case goal_map:
        index_load();
        do_mapping();
        break;
    case goal_unknown:
    default:
        ERROR_EXIT(23, "goal has not been specified!");
    }

    index_free();
    genome_free();
    free_args();

    return 0;
}
