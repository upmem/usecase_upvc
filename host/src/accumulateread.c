/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _GNU_SOURCE

#include "accumulateread.h"
#include "common.h"
#include "index.h"
#include "upvc.h"
#include "genome.h"

#include <assert.h>
#include <dpu_backend.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>
#include <unistd.h>
#include <limits.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))

static FILE **result_file;
static acc_results_t *results_buffers[NB_DISPATCH_AND_ACC_BUFFER];
#define RESULTS_BUFFERS(pass_id) results_buffers[(pass_id) % NB_DISPATCH_AND_ACC_BUFFER]

#define BUCKET_SIZE (16)
#define NB_BUCKET (1 << BUCKET_SIZE)
#define BUCKET_MASK (NB_BUCKET - 1)
typedef struct bucket_elem {
    dpu_result_out_t *elem;
    struct bucket_elem *next;
} bucket_elem_t;


#define ACCUMULATE_THREADS (8)
#define ACCUMULATE_THREADS_SLAVE (ACCUMULATE_THREADS - 1)
static acc_results_t *acc_res;

acc_results_t accumulate_get_result(unsigned int pass_id)
{
    if (result_file[pass_id] == NULL) {
        static const dpu_result_out_t dummy_res = { .num = -1 };
        char result_filename[512];
        sprintf(result_filename, "result_%u.bin", pass_id);

        result_file[pass_id] = fopen(result_filename, "w+");
        assert(result_file[pass_id] != NULL);
        assert(unlink(result_filename) == 0);
        fwrite(&dummy_res, sizeof(dummy_res), 1, result_file[pass_id]);
    }

    fseek(result_file[pass_id], 0, SEEK_END);
    size_t size = ftell(result_file[pass_id]);
    rewind(result_file[pass_id]);

    dpu_result_out_t *results = (dpu_result_out_t *)malloc(size);
    assert(results != NULL);
    size_t size_read = fread(results, size, 1, result_file[pass_id]);
    assert(size_read == 1);

    return (acc_results_t) { .nb_res = (size / sizeof(dpu_result_out_t)) - 1, .results = results };
}

#define NB_RANK_MAX (40)
#define NB_CI (8)
#define NB_DPU_PER_CI (8)

static uint32_t dpu_error_res[NB_RANK_MAX][NB_CI][NB_DPU_PER_CI] = { 0 };
static bool error_in_run = false;

void accumulate_summary_dpu_disabled() {
    printf("\n\n######################################\n"
           "### 3/3 OK! stopping now :) ##########\n\n");

    printf("DISABLED DPUs:\n");
    for (uint32_t each_rank = 0; each_rank < NB_RANK_MAX; each_rank++) {
        for (uint32_t each_ci = 0; each_ci < NB_CI; each_ci++) {
            for (uint32_t each_dpu = 0; each_dpu < NB_DPU_PER_CI; each_dpu++) {
                if (dpu_error_res[each_rank][each_ci][each_dpu] != 0) {
                    printf("\t%u.%u.%u\n", each_rank, each_ci, each_dpu);
                }
            }
        }
    }
    printf("\n######################################\n");
}
static void disable_dpu(uint32_t rank, uint32_t ci, uint32_t dpu, FILE *f_disabled_dpus) {
    char *command;
    printf("Disabling dpu %u.%u.%u (%u errors) ...", rank, ci, dpu, dpu_error_res[rank][ci][dpu]);
    assert(asprintf(&command, "upmem-dimm-configure.py --rank /dev/dpu_rank%u --disable-dpu %u.%u", rank, ci, dpu) > 0);
    int ret = system(command);
    if (ret == 0) {
        printf("OK!\n");
        dpu_error_res[rank][ci][dpu] = UINT_MAX;
        fprintf(f_disabled_dpus, "%u.%u.%u disabled\n", rank, ci, dpu);
    } else {
        printf("ERROR: %u\n", ret);
    }
    free(command);
}

static void print_time(FILE *f_disabled_dpus, int round)
{
    time_t timer;
    char time_buf[26];
    struct tm *tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    fprintf(f_disabled_dpus, "%s run %u\n", time_buf, round);
}


void accumulate_disable_dpus()
{
    FILE *f_disabled_dpus = fopen("disabled_dpus.txt", "a");
    assert(f_disabled_dpus != NULL);
    static int round = 1;
    print_time(f_disabled_dpus, round++);
    printf("### DISABLING DPUS ###########################\n");
    for (uint32_t each_rank = 0; each_rank < NB_RANK_MAX; each_rank++) {
        for (uint32_t each_ci = 0; each_ci < NB_CI; each_ci++) {
            for (uint32_t each_dpu = 0; each_dpu < NB_DPU_PER_CI; each_dpu++) {
                if (dpu_error_res[each_rank][each_ci][each_dpu] != 0 && dpu_error_res[each_rank][each_ci][each_dpu] != UINT_MAX) {
                    disable_dpu(each_rank, each_ci, each_dpu, f_disabled_dpus);
                }
            }
        }
    }
    printf("##############################################\n");
    fclose(f_disabled_dpus);
}

#define NB_CHAR_PER_LINE 50
void accumulate_read(unsigned int pass_id, uint32_t max_nb_pass)
{
    acc_res = RESULTS_BUFFERS(pass_id);
    genome_t *genome = genome_get();
    uint32_t nb_char_max_nb_pass = snprintf(NULL, 0, "%u", max_nb_pass);
    uint32_t done = NB_CHAR_PER_LINE * (float)(pass_id+1)/max_nb_pass;
    static char progress_bar[NB_CHAR_PER_LINE + 1];
    uint32_t each_char = 0;
    for (; each_char < done; each_char++) {
        progress_bar[each_char] = '#';
    }
    for (; each_char < NB_CHAR_PER_LINE; each_char++) {
        progress_bar[each_char] = '-';
    }
    progress_bar[each_char] = '\0';
    printf("\rPASS %*u [%s]", nb_char_max_nb_pass, pass_id, progress_bar);
    fflush(stdout);

    // compute the total number of resultat for all DPUs
    for (unsigned int numdpu = 0; numdpu < nb_dpus_per_run; numdpu++) {
        uint32_t rank, ci, dpu;
        get_dpu_info(numdpu, &rank, &ci, &dpu);
        assert(rank < NB_RANK_MAX && ci < NB_CI && dpu < NB_DPU_PER_CI);
        if (acc_res[numdpu].results[acc_res[numdpu].nb_res].num != -1) {
            uint32_t nb_error = dpu_error_res[rank][ci][dpu]++;
            error_in_run = true;
            if (nb_error == 0) {
                printf("\n!! DPU %u.%u.%u wrong end-mark !!\n", rank, ci, dpu);
            }
        }
        for (unsigned int each_res = 0; each_res < acc_res[numdpu].nb_res; each_res++) {
            uint32_t seq_nr = acc_res[numdpu].results[each_res].coord.seq_nr;
            if ((seq_nr >= genome->nb_seq) || (acc_res[numdpu].results[each_res].coord.seed_nr >= genome->len_seq[seq_nr])){
                uint32_t nb_error = dpu_error_res[rank][ci][dpu]++;
                error_in_run = true;
                if (nb_error == 0) {
                    printf("\n!! DPU %u.%u.%u wrong result (%u/%u, %u/%lu) !!\n", rank, ci, dpu, seq_nr, genome->nb_seq,
                        acc_res[numdpu].results[each_res].coord.seed_nr, genome->len_seq[seq_nr]);
                }
            }
        }
    }
}

uint32_t accumulate_valid_run()
{
    static uint32_t nb_run_valid = 0;
    if (!error_in_run) {
        ++nb_run_valid;
    } else {
        nb_run_valid = 0;
    }
    error_in_run = false;

    return nb_run_valid;
}

void accumulate_free()
{
    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
            free(results_buffers[each_pass][each_dpu].results);
        }
        free(results_buffers[each_pass]);
    }
}

void accumulate_init()
{
    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        results_buffers[each_pass] = (acc_results_t *)malloc(sizeof(acc_results_t) * nb_dpus_per_run);
        assert(results_buffers[each_pass] != NULL);
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
            results_buffers[each_pass][each_dpu].results = (dpu_result_out_t *)malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
            assert(results_buffers[each_pass][each_dpu].results != NULL);
        }
    }
}

acc_results_t *accumulate_get_buffer(unsigned int dpu_id, unsigned int pass_id) { return &(RESULTS_BUFFERS(pass_id)[dpu_id]); }
