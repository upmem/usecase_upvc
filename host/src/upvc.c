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
#include "debug.h"
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
#include "profiling.h"

#include "backends_functions.h"

unsigned int nb_dpus_per_run;

static backends_functions_t backends_functions;
static unsigned int round;
static FILE *fipe1, *fipe2, *fope1, *fope2;
static sem_t getreads_to_dispatch_sem, dispatch_to_exec_sem, exec_to_dispatch_sem, exec_to_acc_sem, acc_to_exec_sem,
    accprocess_to_getreads_sem, acc_to_process_sem;

#define LAST_RUN(dpu_offset) (((dpu_offset) + nb_dpus_per_run) >= index_get_nb_dpu())
#define FOREACH_RUN(dpu_offset) for (unsigned int dpu_offset = 0; dpu_offset < index_get_nb_dpu(); dpu_offset += nb_dpus_per_run)
#define FOREACH_PASS(each_pass) for (unsigned int each_pass = 0; get_reads_in_buffer(each_pass) != 0; each_pass++)
#define FOR(loop) for (int _i = 0; _i < (loop); _i++)

void *thread_get_reads(__attribute__((unused)) void *arg)
{
    STAT_RECORD_START(STAT_THREAD_GET_READS);
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_READS_BUFFER) { sem_post(&accprocess_to_getreads_sem); }

        unsigned int each_pass = 0;
        do {
            sem_wait(&accprocess_to_getreads_sem);
            get_reads(fipe1, fipe2, each_pass);
            sem_post(&getreads_to_dispatch_sem);
        } while (get_reads_in_buffer(each_pass++) != 0);

        fseek(fipe1, 0, SEEK_SET);
        fseek(fipe2, 0, SEEK_SET);

        FOR(NB_READS_BUFFER) { sem_wait(&accprocess_to_getreads_sem); }
    }
    STAT_RECORD_LAST_STEP(STAT_THREAD_GET_READS, 0);

    return NULL;
}

void exec_dpus()
{
    STAT_RECORD_START(STAT_EXEC_DPUS);
    FOREACH_RUN(dpu_offset)
    {
        backends_functions.load_mram(dpu_offset, 0);

        STAT_RECORD_STEP(STAT_EXEC_DPUS, 0);

        sem_wait(&dispatch_to_exec_sem);

        STAT_RECORD_STEP(STAT_EXEC_DPUS, 1);

        FOREACH_PASS(each_pass)
        {
            backends_functions.run_dpu(
                dpu_offset, each_pass, &exec_to_dispatch_sem, &acc_to_exec_sem, &exec_to_acc_sem, &dispatch_to_exec_sem);
        }

        STAT_RECORD_STEP(STAT_EXEC_DPUS, 2);

        backends_functions.wait_dpu();

        STAT_RECORD_STEP(STAT_EXEC_DPUS, 3);

        sem_post(&exec_to_acc_sem);
    }
    STAT_RECORD_LAST_STEP(STAT_EXEC_DPUS, 4);
}

void *thread_dispatch(__attribute__((unused)) void *arg)
{
    STAT_RECORD_START(STAT_THREAD_DISPATCH);
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            sem_post(&exec_to_dispatch_sem);
        }

        sem_wait(&getreads_to_dispatch_sem);
        FOREACH_PASS(each_pass)
        {
            sem_wait(&exec_to_dispatch_sem);
            dispatch_read(each_pass);
            sem_post(&dispatch_to_exec_sem);
            sem_wait(&getreads_to_dispatch_sem);
        }

        sem_post(&dispatch_to_exec_sem);

        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            sem_wait(&exec_to_dispatch_sem);
        }
    }
    STAT_RECORD_LAST_STEP(STAT_THREAD_DISPATCH, 0);
    return NULL;
}

void *thread_acc(__attribute__((unused)) void *arg)
{
    STAT_RECORD_START(STAT_THREAD_ACC);
    FOREACH_RUN(dpu_offset)
    {
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            sem_post(&acc_to_exec_sem);
        }

        sem_wait(&exec_to_acc_sem);

        FOREACH_PASS(each_pass)
        {
            accumulate_read(each_pass, dpu_offset);

            sem_post(&acc_to_exec_sem);
            if (LAST_RUN(dpu_offset)) {
                sem_post(&acc_to_process_sem);
            } else {
                sem_post(&accprocess_to_getreads_sem);
            }
            sem_wait(&exec_to_acc_sem);
        }
        if (LAST_RUN(dpu_offset)) {
            sem_post(&acc_to_process_sem);
        } else {
            sem_post(&accprocess_to_getreads_sem);
        }
        FOR(NB_DISPATCH_AND_ACC_BUFFER)
        {
            sem_wait(&acc_to_exec_sem);
        }
    }
    STAT_RECORD_LAST_STEP(STAT_THREAD_ACC, 0);
    return NULL;
}

void *thread_process(__attribute__((unused)) void *arg)
{
    struct timespec start_time, current_time;
    STAT_RECORD_START(STAT_THREAD_PROCESS);
    sem_wait(&acc_to_process_sem);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
    STAT_RECORD_STEP(STAT_THREAD_PROCESS, 0);

    FOREACH_PASS(each_pass)
    {
        process_read(fope1, fope2, round, each_pass);
        clock_gettime(CLOCK_MONOTONIC_RAW, &current_time);
        int spent = current_time.tv_sec-start_time.tv_sec;
        printf("time spent: %02dh%02dm%02ds (%ds)\n", spent/3600, spent/60%60, spent%60, spent);
        sem_post(&accprocess_to_getreads_sem);
        sem_wait(&acc_to_process_sem);
    }

    STAT_RECORD_STEP(STAT_THREAD_PROCESS, 1);
    sem_post(&accprocess_to_getreads_sem);
    STAT_RECORD_LAST_STEP(STAT_THREAD_PROCESS, 2);

    return NULL;
}

static void exec_round()
{
    STAT_RECORD_START(STAT_EXEC_ROUND);
    char filename[FILENAME_MAX];
    char *input_prefix = get_input_path();
    static unsigned int max_nb_pass;

    if (round == 0) {
        size_t read_size1, read_size2, nb_read1, nb_read2;

        sprintf(filename, "%s_PE1.fastq", input_prefix);
        fipe1 = fopen(filename, "r");
        CHECK_FILE(fipe1, filename);
        assert(get_input_info(fipe1, &read_size1, &nb_read1) == 0);
        printf("read_size1 %zu SIZE_READ %u\n", read_size1, SIZE_READ);
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

    STAT_RECORD_STEP(STAT_EXEC_ROUND, 0);
    accumulate_init(max_nb_pass);
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 1);

    pthread_t tid_get_reads;
    pthread_t tid_dispatch;
    pthread_t tid_acc;
    pthread_t tid_process;

    int ret;

    // INIT
    ret = sem_init(&getreads_to_dispatch_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&dispatch_to_exec_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&exec_to_dispatch_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&exec_to_acc_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&acc_to_exec_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&acc_to_process_sem, 0, 0);
    assert(ret == 0);
    ret = sem_init(&accprocess_to_getreads_sem, 0, 0);
    assert(ret == 0);
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 2);

    // CREATE
    ret = pthread_create(&tid_get_reads, NULL, thread_get_reads, NULL);
    assert(ret == 0);
    ret = pthread_create(&tid_dispatch, NULL, thread_dispatch, NULL);
    assert(ret == 0);
    ret = pthread_create(&tid_acc, NULL, thread_acc, NULL);
    assert(ret == 0);
    ret = pthread_create(&tid_process, NULL, thread_process, NULL);
    assert(ret == 0);
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 3);

    // EXECUTE
    exec_dpus();
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 4);

    // JOIN
    ret = pthread_join(tid_get_reads, NULL);
    assert(ret == 0);
    ret = pthread_join(tid_dispatch, NULL);
    assert(ret == 0);
    ret = pthread_join(tid_acc, NULL);
    assert(ret == 0);
    ret = pthread_join(tid_process, NULL);
    assert(ret == 0);
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 5);

    // DESTROY
    ret = sem_destroy(&getreads_to_dispatch_sem);
    assert(ret == 0);
    ret = sem_destroy(&dispatch_to_exec_sem);
    assert(ret == 0);
    ret = sem_destroy(&exec_to_dispatch_sem);
    assert(ret == 0);
    ret = sem_destroy(&acc_to_exec_sem);
    assert(ret == 0);
    ret = sem_destroy(&exec_to_acc_sem);
    assert(ret == 0);
    ret = sem_destroy(&acc_to_process_sem);
    assert(ret == 0);
    ret = sem_destroy(&accprocess_to_getreads_sem);
    assert(ret == 0);
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 6);

    accumulate_free();
    STAT_RECORD_STEP(STAT_EXEC_ROUND, 7);

    fclose(fipe1);
    fclose(fipe2);
    STAT_RECORD_LAST_STEP(STAT_EXEC_ROUND, 8);
}

static void do_mapping()
{
    STAT_RECORD_START(STAT_DO_MAPPING);
    backends_functions.init_backend(&nb_dpus_per_run);
    variant_tree_init();
    dispatch_init();
    process_read_init();
    STAT_RECORD_STEP(STAT_DO_MAPPING, 0);

    for (round = 0; round < NB_ROUND; round++) {
        printf("#################\n"
               "starting round %u\n"
               "#################\n",
            round);
        exec_round();
    }
    STAT_RECORD_STEP(STAT_DO_MAPPING, 1);
    create_vcf();
    STAT_RECORD_STEP(STAT_DO_MAPPING, 2);

    process_read_free();
    dispatch_free();
    variant_tree_free();
    backends_functions.free_backend();
    STAT_RECORD_LAST_STEP(STAT_DO_MAPPING, 3);
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
    memset(profiling, 0, sizeof(profiling));

    printf("%s\n", VERSION);
    print_time();
    for (int i=0; i<15; i++)
            profiling[i] = (struct time_stat_t){0};

    printf("Information:\n");
    printf("\tread size: %d\n", SIZE_READ);
    printf("\tseed size: %u\n", SIZE_SEED);
    struct timespec start_time, start_process_time, stop_time, stop_process_time;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_process_time);

    if (get_simulation_mode()) {
        backends_functions.init_backend = init_backend_simulation;
        backends_functions.free_backend = free_backend_simulation;
        backends_functions.run_dpu = run_dpu_simulation;
        backends_functions.load_mram = load_mram_simulation;
        backends_functions.wait_dpu = wait_dpu_simulation;
    } else {
        backends_functions.init_backend = init_backend_dpu;
        backends_functions.free_backend = free_backend_dpu;
        backends_functions.run_dpu = run_on_dpu;
        backends_functions.load_mram = load_mram_dpu;
        backends_functions.wait_dpu = wait_dpu_dpu;
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

    clock_gettime(CLOCK_MONOTONIC_RAW, &stop_time);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_process_time);

    double time = (float)((stop_time.tv_sec - start_time.tv_sec) * 1e9 + stop_time.tv_nsec - start_time.tv_nsec) / (1e9);
    double process_time = (float)((stop_process_time.tv_sec - start_process_time.tv_sec) * 1e9 + stop_process_time.tv_nsec
                              - start_process_time.tv_nsec)
        / (1e9);
    printf("time: %f\n"
           "process time: %f\n"
           "ratio: %f\n",
        time, process_time, process_time / time);

    printf("\nprofiling:\n\n");
    PRINT_FUNCTION_STAT(STAT_DO_MAPPING);
    PRINT_FUNCTION_STAT(STAT_THREAD_PROCESS);
    PRINT_FUNCTION_STAT(STAT_THREAD_ACC);
    PRINT_FUNCTION_STAT(STAT_THREAD_DISPATCH);
    PRINT_FUNCTION_STAT(STAT_THREAD_GET_READS);
    PRINT_FUNCTION_STAT(STAT_EXEC_DPUS);
    PRINT_FUNCTION_STAT(STAT_EXEC_ROUND);
    PRINT_FUNCTION_STAT(STAT_PROCESS_READ);
    PRINT_FUNCTION_STAT(STAT_DO_PROCESS_READ);
    PRINT_FUNCTION_STAT(STAT_ADD_TO_NON_MAPPED_READ);
    PRINT_FUNCTION_STAT(STAT_UPDATE_FREQUENCY_TABLE);
    PRINT_FUNCTION_STAT(STAT_GET_READ_UPDATE_POSITIONS);
    PRINT_FUNCTION_STAT(STAT_CODE_ALIGNMENT);
    PRINT_FUNCTION_STAT(STAT_SUB_ONLY_PATH);
    PRINT_FUNCTION_STAT(STAT_DPD);

    index_free();
    genome_free();
    free_args();

    return 0;
}
