/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <semaphore.h>

#include "parse_args.h"
#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "processread.h"
#include "vartree.h"
#include "vcf.h"
#include "genome.h"
#include "upvc_dpu.h"
#include "upvc.h"
#include "mram_dpu.h"
#include "backends_functions.h"
#include "simu_backend.h"
#include "dpu_backend.h"

#define MAX_NB_PASS (4096)
#define NB_READS_BUFFER (8)

static volatile unsigned int global_thread_process_each_pass_mod_start;

typedef struct {
        sem_t *dispatch_free_sem;
        sem_t *acc_process_wait_sem;
        int8_t **reads_buffer;
        int *nb_read;
        FILE *fipe1;
        FILE *fipe2;
        times_ctx_t *times_ctx;
        reads_info_t *reads_info;
} thread_get_reads_arg_t;

void *thread_get_reads(void *arg)
{
        thread_get_reads_arg_t *args = (thread_get_reads_arg_t *)arg;
        int8_t **reads_buffer = args->reads_buffer;
        int *nb_read = args->nb_read;
        FILE *fipe1 = args->fipe1;
        FILE *fipe2 = args->fipe2;
        sem_t *dispatch_free_sem = args->dispatch_free_sem;
        sem_t *acc_process_wait_sem = args->acc_process_wait_sem;
        times_ctx_t *times_ctx = args->times_ctx;
        reads_info_t *reads_info = args->reads_info;

        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        unsigned int each_pass_mod = 0;

        for (unsigned int dpu_offset = 0; dpu_offset < nb_dpu; dpu_offset += nb_dpus_per_run) {
                unsigned int each_pass = 0;

                sem_wait(acc_process_wait_sem);

                PRINT_TIME_GET_READS(times_ctx, each_pass);
                reads_buffer[each_pass_mod] = (int8_t *) malloc(sizeof(int8_t) * MAX_READS_BUFFER * reads_info->size_read);
                nb_read[each_pass_mod] = get_reads(fipe1, fipe2, reads_buffer[each_pass_mod], times_ctx, reads_info);
                PRINT_TIME_GET_READS(times_ctx, each_pass);

                sem_post(dispatch_free_sem);

                while (nb_read[each_pass_mod] != 0) {
                        if (DEBUG_PASS != -1) {
                                if ((int)each_pass == DEBUG_PASS) {
                                        each_pass++;
                                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
                                        nb_read[each_pass_mod] = 0;
                                        sem_post(dispatch_free_sem);
                                        break;
                                } else {
                                        continue;
                                }
                        }

                        each_pass++;
                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
                        assert(each_pass < MAX_NB_PASS);

                        sem_wait(acc_process_wait_sem);

                        PRINT_TIME_GET_READS(times_ctx, each_pass);
                        reads_buffer[each_pass_mod] = (int8_t *) malloc(sizeof(int8_t) * MAX_READS_BUFFER * reads_info->size_read);
                        nb_read[each_pass_mod] = get_reads(fipe1, fipe2, reads_buffer[each_pass_mod], times_ctx, reads_info);
                        PRINT_TIME_GET_READS(times_ctx, each_pass);

                        sem_post(dispatch_free_sem);
                }
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
        }

        return NULL;
}

typedef struct {
        sem_t *dispatch_wait_sem;
        sem_t *dispatch_free_sem;
        sem_t *acc_wait_sem;
        sem_t *acc_free_sem;
        unsigned int rank_id;
        int *nb_read;
        devices_t *devices;
        dispatch_request_t *dispatch_requests;
        backends_functions_t *backends_functions;
        times_ctx_t *times_ctx;
        reads_info_t *reads_info;
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
        volatile int *nb_read = args->nb_read;
        devices_t *devices = args->devices;
        dispatch_request_t *dispatch_requests = args->dispatch_requests;
        backends_functions_t *backends_functions = args->backends_functions;
        times_ctx_t *times_ctx = args->times_ctx;
        reads_info_t *reads_info = args->reads_info;

        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        unsigned int each_pass_mod = 0;

        for (unsigned int dpu_offset = 0; dpu_offset < nb_dpu; dpu_offset += nb_dpus_per_run) {
                printf("[%u-%u] loading mram\n", dpu_offset, rank_id);
                PRINT_TIME_WRITE_MRAM(times_ctx, 0, rank_id);
                t1 = my_clock();
                backends_functions->load_mram(dpu_offset, rank_id, devices, reads_info, times_ctx);
                t2 = my_clock();
                PRINT_TIME_WRITE_MRAM(times_ctx, 0, rank_id);
                printf("[%u-%u]  time %.2lf sec\n", dpu_offset, rank_id, t2 - t1);

                unsigned int each_pass = 0;
                if (DEBUG_PASS != -1) {
                        each_pass = DEBUG_PASS;
                }

                sem_wait(dispatch_wait_sem);
                while (nb_read[each_pass_mod] != 0) {

                        printf("[%u-%u] running pass\n", each_pass, rank_id);
                        PRINT_TIME_MAP_READ(times_ctx, each_pass, rank_id);
                        t1 = my_clock();
                        backends_functions->run_dpu(dispatch_requests,
                                                    devices,
                                                    dpu_offset,
                                                    rank_id,
                                                    each_pass,
                                                    dispatch_free_sem,
                                                    acc_wait_sem,
                                                    times_ctx,
                                                    reads_info);
                        t2 = my_clock();
                        PRINT_TIME_MAP_READ(times_ctx, each_pass, rank_id);
                        printf("[%u-%u]  time %.2lf sec\n", each_pass, rank_id, t2 - t1);

                        sem_post(acc_free_sem);
                        each_pass++;
                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
                        sem_wait(dispatch_wait_sem);
                }
                sem_post(acc_free_sem);
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
        }
        return NULL;
}

typedef struct {
        sem_t *get_reads_wait_sem;
        sem_t *exec_rank_wait_sem;
        sem_t *exec_rank_free_sem;
        unsigned int nb_rank;
        int *nb_read;
        int8_t **reads_buffer;
        index_seed_t **index_seed;
        dispatch_request_t *dispatch_requests;
        backends_functions_t *backends_functions;
        times_ctx_t *times_ctx;
        reads_info_t *reads_info;
} thread_dispatch_arg_t;

void *thread_dispatch(void *arg)
{
        thread_dispatch_arg_t *args = (thread_dispatch_arg_t *)arg;
        sem_t *get_reads_wait_sem = args->get_reads_wait_sem;
        sem_t *exec_rank_wait_sem = args->exec_rank_wait_sem;
        sem_t *exec_rank_free_sem = args->exec_rank_free_sem;
        unsigned int nb_rank = args->nb_rank;
        volatile int *nb_read = args->nb_read;
        int8_t **reads_buffer = args->reads_buffer;
        index_seed_t **index_seed = args->index_seed;
        dispatch_request_t *dispatch_requests = args->dispatch_requests;
        backends_functions_t *backends_functions = args->backends_functions;
        times_ctx_t *times_ctx = args->times_ctx;
        reads_info_t *reads_info = args->reads_info;

        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        unsigned int each_pass_mod = 0;

        for (unsigned int dpu_offset = 0; dpu_offset < nb_dpu; dpu_offset += nb_dpus_per_run) {
                unsigned int each_pass = 0;
                if (DEBUG_PASS != -1) {
                        each_pass = DEBUG_PASS;
                }

                sem_wait(get_reads_wait_sem);

                while (nb_read[each_pass_mod] != 0) {
                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_wait(&exec_rank_wait_sem[each_rank]);
                        }

                        PRINT_TIME_DISPATCH(times_ctx, each_pass);
                        dispatch_read(index_seed,
                                      reads_buffer[each_pass_mod],
                                      nb_read[each_pass_mod],
                                      dispatch_requests,
                                      times_ctx,
                                      reads_info,
                                      backends_functions);
                        PRINT_TIME_DISPATCH(times_ctx, each_pass);

                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_post(&exec_rank_free_sem[each_rank]);
                        }

                        each_pass++;
                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
                        sem_wait(get_reads_wait_sem);
                }
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        sem_post(&exec_rank_free_sem[each_rank]);
                }
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
        }
        return NULL;
}

typedef struct {
        sem_t *exec_rank_wait_sem;
        sem_t *exec_rank_free_sem;
        sem_t *process_free_sem;
        sem_t *get_reads_free_sem;
        unsigned int nb_rank;
        int *nb_read;
        int8_t **reads_buffer;
        unsigned int *result_tab_nb_read;
        dpu_result_out_t **result_tab;
        times_ctx_t *times_ctx;
} thread_acc_arg_t;

void *thread_acc(void *arg)
{
        thread_acc_arg_t *args = (thread_acc_arg_t*)arg;
        sem_t *exec_rank_wait_sem = args->exec_rank_wait_sem;
        sem_t *exec_rank_free_sem = args->exec_rank_free_sem;
        sem_t *process_free_sem = args->process_free_sem;
        sem_t *get_reads_free_sem = args->get_reads_free_sem;
        unsigned int nb_rank = args->nb_rank;
        volatile int *nb_read = args->nb_read;
        int8_t **reads_buffer = args->reads_buffer;
        unsigned int *result_tab_nb_read = args->result_tab_nb_read;
        dpu_result_out_t **result_tab = args->result_tab;
        times_ctx_t *times_ctx = args->times_ctx;

        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        unsigned int dpu_offset = 0;
        unsigned int each_pass_mod = 0;

        for (; dpu_offset < (nb_dpu - nb_dpus_per_run); dpu_offset += nb_dpus_per_run) {
                unsigned int each_pass = 0;
                if (DEBUG_PASS != -1) {
                        each_pass = DEBUG_PASS;
                }

                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        sem_wait(&exec_rank_wait_sem[each_rank]);
                }
                while (nb_read[each_pass] != 0) {

                        PRINT_TIME_ACC_READ(times_ctx, each_pass);
                        accumulate_read(&result_tab[each_pass], &result_tab_nb_read[each_pass], dpu_offset, times_ctx);
                        PRINT_TIME_ACC_READ(times_ctx, each_pass);

                        free(reads_buffer[each_pass_mod]);
                        sem_post(get_reads_free_sem);
                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_post(&exec_rank_free_sem[each_rank]);
                        }
                        each_pass++;
                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;

                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_wait(&exec_rank_wait_sem[each_rank]);
                        }
                }
                sem_post(get_reads_free_sem);
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
        }
        {
                unsigned int each_pass = 0;
                global_thread_process_each_pass_mod_start = each_pass_mod;
                if (DEBUG_PASS != -1) {
                        each_pass = DEBUG_PASS;
                }

                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        sem_wait(&exec_rank_wait_sem[each_rank]);
                }
                while (nb_read[each_pass_mod] != 0) {

                        PRINT_TIME_ACC_READ(times_ctx, each_pass);
                        accumulate_read(&result_tab[each_pass], &result_tab_nb_read[each_pass], dpu_offset, times_ctx);
                        PRINT_TIME_ACC_READ(times_ctx, each_pass);

                        sem_post(process_free_sem);
                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_post(&exec_rank_free_sem[each_rank]);
                        }
                        each_pass++;
                        each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;

                        for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                                sem_wait(&exec_rank_wait_sem[each_rank]);
                        }
                }
                sem_post(process_free_sem);
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
        }
        return NULL;
}

typedef struct {
        sem_t *acc_wait_sem;
        sem_t *get_reads_free_sem;
        int round;
        int *nb_read;
        genome_t *ref_genome;
        int8_t **reads_buffer;
        variant_tree_t **variant_list;
        int *substitution_list;
        int8_t *mapping_coverage;
        unsigned int *result_tab_nb_read;
        dpu_result_out_t **result_tab;
        FILE *fope1;
        FILE *fope2;
        times_ctx_t *times_ctx;
        reads_info_t *reads_info;
} thread_process_arg_t;

void *thread_process(void *arg)
{
        thread_process_arg_t *args = (thread_process_arg_t *)arg;
        sem_t *acc_wait_sem = args->acc_wait_sem;
        sem_t *get_reads_free_sem = args->get_reads_free_sem;
        int round = args->round;
        volatile int *nb_read = args->nb_read;
        genome_t *ref_genome = args->ref_genome;
        int8_t **reads_buffer = args->reads_buffer;
        variant_tree_t **variant_list = args->variant_list;
        int *substitution_list = args->substitution_list;
        int8_t *mapping_coverage = args->mapping_coverage;
        unsigned int *result_tab_nb_read = args->result_tab_nb_read;
        dpu_result_out_t **result_tab = args->result_tab;
        FILE *fope1 = args->fope1;
        FILE *fope2 = args->fope2;
        times_ctx_t *times_ctx = args->times_ctx;
        reads_info_t *reads_info = args->reads_info;

        unsigned int each_pass = 0, each_pass_mod;
        if (DEBUG_PASS != -1) {
                each_pass = DEBUG_PASS;
        }

        sem_wait(acc_wait_sem);
        each_pass_mod = global_thread_process_each_pass_mod_start;
        while (nb_read[each_pass_mod] != 0) {
                PRINT_TIME_PROCESS_READ(times_ctx, each_pass);
                process_read(ref_genome,
                             reads_buffer[each_pass_mod],
                             variant_list,
                             substitution_list,
                             mapping_coverage,
                             result_tab[each_pass],
                             result_tab_nb_read[each_pass],
                             fope1,
                             fope2,
                             round,
                             times_ctx,
                             reads_info);
                free(reads_buffer[each_pass_mod]);
                PRINT_TIME_PROCESS_READ(times_ctx, each_pass);
                sem_post(get_reads_free_sem);

                each_pass++;
                each_pass_mod = (each_pass_mod + 1) % NB_READS_BUFFER;
                sem_wait(acc_wait_sem);
        }
        return NULL;
}

static void exec_round(unsigned int round,
                       unsigned int nb_rank,
                       int8_t *mapping_coverage,
                       int *substitution_list,
                       dispatch_request_t *dispatch_requests,
                       index_seed_t **index_seed,
                       char *input_prefix,
                       variant_tree_t **variant_list,
                       genome_t *ref_genome,
                       devices_t *devices,
                       times_ctx_t *times_ctx,
                       reads_info_t *reads_info,
                       backends_functions_t * backends_functions)
{
        char filename[1024];
        FILE *fipe1, *fipe2, *fope1, *fope2;
        int8_t *reads_buffer[NB_READS_BUFFER];
        int nb_read[NB_READS_BUFFER];
        unsigned int result_tab_nb_read[MAX_NB_PASS];
        dpu_result_out_t *result_tab[MAX_NB_PASS];

        memset(result_tab, 0, sizeof(dpu_result_out_t *) * MAX_NB_PASS);
        memset(result_tab_nb_read, 0, sizeof(unsigned int) * MAX_NB_PASS);

        reads_info->delta_neighbour_in_bytes = (SIZE_SEED * round)/4;
        reads_info->size_neighbour_in_32bits_words =
                (reads_info->size_neighbour_in_bytes-reads_info->delta_neighbour_in_bytes) * 4;

        sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round+1);
        fope1 = fopen(filename, "w");
        sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round+1);
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

        {
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

                thread_get_reads_arg_t get_reads_arg =
                        {
                         .dispatch_free_sem = &get_reads_dispatch_sem,
                         .acc_process_wait_sem = &acc_process_get_reads_sem,
                         .reads_buffer = reads_buffer,
                         .nb_read = nb_read,
                         .fipe1 = fipe1,
                         .fipe2 = fipe2,
                         .times_ctx = times_ctx,
                         .reads_info = reads_info,
                        };
                thread_exec_rank_arg_t exec_rank_arg[nb_rank];
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        exec_rank_arg[each_rank].dispatch_wait_sem = &dispatch_exec_rank_sem[each_rank];
                        exec_rank_arg[each_rank].dispatch_free_sem = &exec_rank_dispatch_sem[each_rank];
                        exec_rank_arg[each_rank].acc_wait_sem = &acc_exec_rank_sem[each_rank];
                        exec_rank_arg[each_rank].acc_free_sem = &exec_rank_acc_sem[each_rank];
                        exec_rank_arg[each_rank].rank_id = each_rank;
                        exec_rank_arg[each_rank].nb_read = nb_read;
                        exec_rank_arg[each_rank].devices = devices;
                        exec_rank_arg[each_rank].dispatch_requests = dispatch_requests;
                        exec_rank_arg[each_rank].backends_functions = backends_functions;
                        exec_rank_arg[each_rank].times_ctx = times_ctx;
                        exec_rank_arg[each_rank].reads_info = reads_info;
                }
                thread_dispatch_arg_t dispatch_arg =
                        {
                         .get_reads_wait_sem = &get_reads_dispatch_sem,
                         .exec_rank_wait_sem = exec_rank_dispatch_sem,
                         .exec_rank_free_sem = dispatch_exec_rank_sem,
                         .nb_rank = nb_rank,
                         .nb_read = nb_read,
                         .reads_buffer = reads_buffer,
                         .index_seed = index_seed,
                         .dispatch_requests = dispatch_requests,
                         .backends_functions = backends_functions,
                         .times_ctx = times_ctx,
                         .reads_info = reads_info,
                        };
                thread_acc_arg_t acc_arg =
                        {
                         .exec_rank_wait_sem = exec_rank_acc_sem,
                         .exec_rank_free_sem = acc_exec_rank_sem,
                         .process_free_sem = &acc_process_sem,
                         .get_reads_free_sem = &acc_process_get_reads_sem,
                         .nb_rank = nb_rank,
                         .nb_read = nb_read,
                         .reads_buffer = reads_buffer,
                         .result_tab_nb_read = result_tab_nb_read,
                         .result_tab = result_tab,
                         .times_ctx = times_ctx,
                        };
                thread_process_arg_t process_arg =
                        {
                         .acc_wait_sem = &acc_process_sem,
                         .get_reads_free_sem = &acc_process_get_reads_sem,
                         .round = round,
                         .nb_read = nb_read,
                         .ref_genome = ref_genome,
                         .reads_buffer = reads_buffer,
                         .variant_list = variant_list,
                         .substitution_list = substitution_list,
                         .mapping_coverage = mapping_coverage,
                         .result_tab_nb_read = result_tab_nb_read,
                         .result_tab = result_tab,
                         .fope1 = fope1,
                         .fope2 = fope2,
                         .times_ctx = times_ctx,
                         .reads_info = reads_info,
                        };
                int ret;
                ret = sem_init(&get_reads_dispatch_sem, 0, 0);
                assert(ret == 0);
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
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
                ret = sem_init(&acc_process_get_reads_sem, 0, NB_READS_BUFFER);
                assert(ret == 0);

                ret = pthread_create(&tid_get_reads, NULL, thread_get_reads, (void *)&get_reads_arg);
                assert(ret == 0);
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        ret = pthread_create(&tid_exec_rank[each_rank], NULL, thread_exec_rank, (void *)&exec_rank_arg[each_rank]);
                        assert(ret == 0);
                }
                ret = pthread_create(&tid_dispatch, NULL, thread_dispatch, (void *)&dispatch_arg);
                assert(ret == 0);
                if (DEBUG_DPU == -1) {
                        ret = pthread_create(&tid_acc, NULL, thread_acc, (void *)&acc_arg);
                        assert(ret == 0);
                        ret = pthread_create(&tid_process, NULL, thread_process, (void *)&process_arg);
                        assert(ret == 0);
                }

                ret = pthread_join(tid_get_reads, NULL);
                assert(ret == 0);
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
                        ret = pthread_join(tid_exec_rank[each_rank], NULL);
                        assert(ret == 0);
                }
                ret = pthread_join(tid_dispatch, NULL);
                assert(ret == 0);
                if (DEBUG_DPU == -1) {
                        ret = pthread_join(tid_acc, NULL);
                        assert(ret == 0);
                        ret = pthread_join(tid_process, NULL);
                        assert(ret == 0);
                }

                ret = sem_destroy(&get_reads_dispatch_sem);
                assert(ret == 0);
                for (unsigned int each_rank = 0; each_rank < nb_rank; each_rank++) {
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
        }

        fclose(fipe1);
        fclose(fipe2);
        fclose(fope1);
        fclose(fope2);
}

static void reload_and_verify_mram_images(reads_info_t *reads_info)
{
        index_seed_t **index_seed;
        FILE *seed_file;
        unsigned int nb_dpu = get_nb_dpu();
        malloc_dpu(reads_info, nb_dpu);
        index_seed = load_index_seeds();
        seed_file = fopen(SEED_FILE_LOG, "w");
        print_index_seeds(index_seed, seed_file, reads_info);
        fclose(seed_file);
        printf("Please check %s to verify that the indexing is OK\n", SEED_FILE_LOG);

        free_index(index_seed);
        free_dpu(nb_dpu);
}

static void load_index_save_genome(reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);
        index_seed_t **index_seed = index_genome(ref_genome,
                                                 get_nb_dpu(),
                                                 times_ctx,
                                                 reads_info);
        save_index_seeds(index_seed);

        free_genome(ref_genome);
        free_index(index_seed);
        free_dpu(get_nb_dpu());
}

static void do_mapping(backends_functions_t *backends_functions, reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_rank;
        index_seed_t **index_seed;
        char *input_prefix = get_input_path();
        variant_tree_t *variant_list = NULL;
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);
        devices_t *devices = NULL;
        char filename[1024];

        int8_t *mapping_coverage = (int8_t *) calloc(sizeof(int8_t), ref_genome->fasta_file_size);
        int *substitution_list = (int *) calloc(sizeof(int), ref_genome->fasta_file_size);
        dispatch_request_t *dispatch_requests = dispatch_create(nb_dpu, reads_info);

        backends_functions->init_backend(&nb_rank,
                                         &devices,
                                         get_nb_dpus_per_run(),
                                         get_dpu_binary(),
                                         &index_seed,
                                         nb_dpu,
                                         reads_info);

        sprintf(filename, "%s_time.csv", input_prefix);
        times_ctx->time_file = fopen(filename, "w");
        fprintf(times_ctx->time_file,
                "time, get_reads, dispatch, accumulate_read, process_read, "
                "write_mram, write_reads, compute, read_result, map_read\n");

        for (unsigned int round = 0; round < 3; round++) {
                if (DEBUG_ROUND != -1 && DEBUG_ROUND != round) {
                        continue;
                }

                printf("starting round %u\n", round);
                if (devices != NULL) {
                        fprintf(devices->log_file, "round %i\n", round);
                }
                exec_round(round,
                           nb_rank,
                           mapping_coverage,
                           substitution_list,
                           dispatch_requests,
                           index_seed,
                           input_prefix,
                           &variant_list,
                           ref_genome,
                           devices,
                           times_ctx,
                           reads_info,
                           backends_functions);

                if (DEBUG_PASS != -1) {
                        FILE *res_file;
                        unsigned int first_dpu = 0;
                        unsigned int last_dpu = nb_dpu;
                        sprintf(filename, "%s_%u_res.txt", input_prefix, round);
                        res_file = fopen(filename, "w");
                        if (DEBUG_DPU != -1) {
                                first_dpu = DEBUG_DPU;
                                last_dpu = DEBUG_DPU + 1;
                        }
                        for (unsigned int numdpu = first_dpu; numdpu < last_dpu; numdpu++) {
                                int k = 0;
                                fprintf(res_file, "dpu %i\n", numdpu);
                                while (read_out_num(numdpu, k) != -1) {
                                        fprintf(res_file, "R: %u %u %llu\n",
                                                read_out_num(numdpu, k),
                                                read_out_score(numdpu, k),
                                                (unsigned long long)read_out_coord(numdpu, k).coord);
                                        k++;
                                }
                        }
                        fclose(res_file);
                }
        }

        fclose(times_ctx->time_file);

        backends_functions->free_backend(devices, nb_dpu);

        create_vcf(input_prefix, ref_genome, &variant_list, substitution_list, mapping_coverage, times_ctx);

        free_variant_tree(variant_list);
        free_genome(ref_genome);
        free_index(index_seed);
        dispatch_free(dispatch_requests, nb_dpu);
        free(substitution_list);
        free(mapping_coverage);
}

static void print_time()
{
        time_t timer;
        char time_buf[26];
        struct tm* tm_info;

        time(&timer);
        tm_info = localtime(&timer);

        strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
        printf("upvc started at: %s\n", time_buf);
}

int main(int argc, char *argv[])
{
        reads_info_t reads_info;
        times_ctx_t times_ctx;
        backends_functions_t backends_functions;

        memset(&times_ctx, 0, sizeof(times_ctx_t));
        pthread_mutex_init(&times_ctx.time_file_mutex, NULL);

        validate_args(argc, argv);

        printf("%s\n", VERSION);
        print_time();

        reads_info.size_read = get_read_size(get_input_pe1());
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        printf("Information\n");
        printf(" - read size: %d\n", reads_info.size_read);

        setup_dpus_for_target_type(get_target_type());

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

        switch(get_goal()) {
        case goal_index:
                load_index_save_genome(&reads_info, &times_ctx);
                break;
        case goal_check:
                reload_and_verify_mram_images(&reads_info);
                break;
        case goal_map:
                do_mapping(&backends_functions, &reads_info, &times_ctx);
                break;
        case goal_unknown:
        default:
                ERROR_EXIT(23, "goal has not been specified!");
        }

        pthread_mutex_destroy(&times_ctx.time_file_mutex);
        free_args();

        return 0;
}
