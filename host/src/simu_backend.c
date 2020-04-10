/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>
#include <stdio.h>

#include "dispatch.h"
#include "dpus_mgmt.h"
#include "genome.h"
#include "index.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "simu_backend.h"
#include "upvc.h"
#include "upvc_dpu.h"

#include "common.h"

/* #define PRINT_NODP_ODPD_TIME */
#ifdef PRINT_NODP_ODPD_TIME

#define TIME_VAR                                                                                                                 \
    unsigned int nodp_call = 0;                                                                                                  \
    unsigned int odpd_call = 0;                                                                                                  \
    double tmp_time;                                                                                                             \
    double nodp_time = 0.0;                                                                                                      \
    double odpd_time = 0.0;
#define GET_TIME                                                                                                                 \
    do {                                                                                                                         \
        tmp_time = my_clock();                                                                                                   \
    } while (0)
#define STORE_NODP_TIME                                                                                                          \
    do {                                                                                                                         \
        nodp_time += (my_clock() - tmp_time);                                                                                    \
        nodp_call++;                                                                                                             \
    } while (0)
#define STORE_ODPD_TIME                                                                                                          \
    do {                                                                                                                         \
        odpd_time += (my_clock() - tmp_time);                                                                                    \
        odpd_call++;                                                                                                             \
    } while (0)
#define WRITE_BACK_TIME(dpu_arg)                                                                                                 \
    do {                                                                                                                         \
        dpu_arg->nodp_call = nodp_call;                                                                                          \
        dpu_arg->odpd_call = odpd_call;                                                                                          \
        dpu_arg->nodp_time = nodp_time;                                                                                          \
        dpu_arg->odpd_time = odpd_time;                                                                                          \
    } while (0)
#else

#define TIME_VAR
#define GET_TIME
#define STORE_NODP_TIME
#define STORE_ODPD_TIME
#define WRITE_BACK_TIME(dpu_arg)

#endif

#define MAX_SCORE 40

/**
 * @brief Arguments needed for the thread that will be offload onto a DPU.
 */
typedef struct align_on_dpu_arg {
    int numdpu;
    double nodp_time;
    double odpd_time;
    unsigned int nodp_call;
    unsigned int odpd_call;
    int delta_neighbour;
} align_on_dpu_arg_t;

static int min(int a, int b) { return a < b ? a : b; }

#define PQD_INIT_VAL (99)
static void ODPD_compute(
    int i, int j, int8_t *s1, int8_t *s2, int *Pppj, int Pppjm, int *Qppj, int Qlpj, int *Dppj, int Dppjm, int Dlpj, int Dlpjm)
{
    int d = Dlpjm;
    int QP;
    *Pppj = min(Dppjm + COST_GAPO, Pppjm + COST_GAPE);
    *Qppj = min(Dlpj + COST_GAPO, Qlpj + COST_GAPE);
    QP = min(*Pppj, *Qppj);
    if (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) != ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3)) {
        d += COST_SUB;
    }
    *Dppj = min(d, QP);
}

/**
 * @brief Compute the alignment distance by dynamical programming on the diagonals of the matrix.
 * Stops when score is greater than max_score
 */
static int ODPD(int8_t *s1, int8_t *s2, int max_score, int size_neighbour_in_symbols)
{
    int matrix_size = size_neighbour_in_symbols + 1;
    int D[2][matrix_size];
    int P[2][matrix_size];
    int Q[2][matrix_size];
    int diagonal = (NB_DIAG / 2) + 1;

    for (int j = 0; j <= diagonal; j++) {
        P[0][j] = PQD_INIT_VAL;
        Q[0][j] = PQD_INIT_VAL;
        D[0][j] = j * COST_SUB;
    }
    P[1][0] = PQD_INIT_VAL;
    Q[1][0] = PQD_INIT_VAL;

    for (int i = 1; i < diagonal; i++) {
        int min_score = PQD_INIT_VAL;
        int pp = i % 2;
        int lp = (i - 1) % 2;
        D[pp][0] = i * COST_SUB;
        for (int j = 1; j < i + diagonal; j++) {
            ODPD_compute(
                i, j, s1, s2, &P[pp][j], P[pp][j - 1], &Q[pp][j], Q[lp][j], &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
            if (D[pp][j] < min_score) {
                min_score = D[pp][j];
            }
        }
        Q[pp][i + diagonal] = PQD_INIT_VAL;
        D[pp][i + diagonal] = PQD_INIT_VAL;
        if (min_score > max_score) {
            return min_score;
        }
    }

    for (int i = diagonal; i < matrix_size - diagonal; i++) {
        int min_score = PQD_INIT_VAL;
        int pp = i % 2;
        int lp = (i - 1) % 2;
        P[pp][i - diagonal] = PQD_INIT_VAL;
        D[pp][i - diagonal] = PQD_INIT_VAL;
        for (int j = i + 1 - diagonal; j < i + diagonal; j++) {
            ODPD_compute(
                i, j, s1, s2, &P[pp][j], P[pp][j - 1], &Q[pp][j], Q[lp][j], &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
            if (D[pp][j] < min_score) {
                min_score = D[pp][j];
            }
        }
        Q[pp][i + diagonal] = PQD_INIT_VAL;
        D[pp][i + diagonal] = PQD_INIT_VAL;
        if (min_score > max_score) {
            return min_score;
        }
    }
    int min_score = PQD_INIT_VAL;
    for (int i = matrix_size - diagonal; i < matrix_size; i++) {
        int pp = i % 2;
        int lp = (i - 1) % 2;
        P[pp][i - diagonal] = PQD_INIT_VAL;
        D[pp][i - diagonal] = PQD_INIT_VAL;
        for (int j = i + 1 - diagonal; j < matrix_size; j++) {
            ODPD_compute(
                i, j, s1, s2, &P[pp][j], P[pp][j - 1], &Q[pp][j], Q[lp][j], &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
        }
        if (D[pp][matrix_size - 1] < min_score)
            min_score = D[pp][matrix_size - 1];
    }
    int i = matrix_size - 1;
    int pp = i % 2;
    for (unsigned int j = i + 1 - diagonal; j < (unsigned int)matrix_size; j++) {
        if (D[pp][j] < min_score) {
            min_score = D[pp][j];
        }
    }

    return min_score;
}

static int translation_table[256] = { 0, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20, 20, 30, 30,
    30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20, 30,
    30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 30,
    40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30,
    30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30,
    30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30,
    30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 20,
    30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
    20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40 };

/**
 * @brief Optimized version of ODPD (if no INDELS)
 * If it detects INDELS, return -1. In this case we will need the run the full ODPD.
 */
static int noDP(int8_t *s1, int8_t *s2, int max_score, int delta_neighbour)
{
    int score = 0;
    int size_neighbour = SIZE_NEIGHBOUR_IN_BYTES;
    for (int i = 0; i < size_neighbour - delta_neighbour; i++) {
        int s_xor = ((int)(s1[i] ^ s2[i])) & 0xFF;
        int s_translated = translation_table[s_xor];
        if (s_translated > COST_SUB) {
            int j = i + 1;
            /* INDELS detection */
            if (j < size_neighbour - delta_neighbour - 3) {
                int s1_val = ((int *)(&s1[j]))[0];
                int s2_val = ((int *)(&s2[j]))[0];
                if (((s1_val ^ (s2_val >> 2)) & 0x3FFFFFFF) == 0) {
                    return -1;
                }
                if (((s1_val ^ (s2_val >> 4)) & 0xFFFFFFF) == 0) {
                    return -1;
                }
                if (((s1_val ^ (s2_val >> 6)) & 0x3FFFFFF) == 0) {
                    return -1;
                }
                if (((s1_val ^ (s2_val >> 8)) & 0xFFFFFF) == 0) {
                    return -1;
                }
                if (((s2_val ^ (s1_val >> 2)) & 0x3FFFFFFF) == 0) {
                    return -1;
                }
                if (((s2_val ^ (s1_val >> 4)) & 0xFFFFFFF) == 0) {
                    return -1;
                }
                if (((s2_val ^ (s1_val >> 6)) & 0x3FFFFFF) == 0) {
                    return -1;
                }
                if (((s2_val ^ (s1_val >> 8)) & 0xFFFFFF) == 0) {
                    return -1;
                }
            }
        }
        score += s_translated;
        if (score > max_score)
            break;
    }
    return score;
}

static void *align_on_dpu(void *arg)
{
    int nb_map = 0;
    align_on_dpu_arg_t *dpu_arg = (align_on_dpu_arg_t *)arg;
    int numdpu = dpu_arg->numdpu;
    int delta_neighbour = dpu_arg->delta_neighbour;
    int size_neighbour = SIZE_NEIGHBOUR_IN_BYTES;
    int size_neighbour_in_symbols = SIZE_IN_SYMBOLS(delta_neighbour);
    int aligned_nbr_size = ALIGN_DPU(size_neighbour);
    int size_nbr_coord = aligned_nbr_size + sizeof(dpu_result_coord_t);
    mem_dpu_t *M = get_mem_dpu(numdpu);
    dpu_result_out_t *M_res = get_mem_dpu_res(numdpu);
    int current_read = 0;
    TIME_VAR;

    M_res[nb_map].num = -1;
    while (M->num[current_read] != -1) {
        int offset = M->offset[current_read]; /* address of the first neighbour */
        int min = MAX_SCORE;
        int nb_map_start = nb_map;
        for (int nb_neighbour = 0; nb_neighbour < M->count[current_read]; nb_neighbour++) {
            dpu_result_coord_t curr_coord
                = *((dpu_result_coord_t *)(&M->neighbour_idx[(offset + nb_neighbour) * size_nbr_coord]));
            int8_t *curr_nbr = &M->neighbour_idx[(offset + nb_neighbour) * size_nbr_coord + sizeof(dpu_result_coord_t)];
            int8_t *curr_read = &M->neighbour_read[current_read * size_neighbour];

            GET_TIME;
            int score_noDP = noDP(curr_read, curr_nbr, min, delta_neighbour);
            STORE_NODP_TIME;

            int score = score_noDP;
            if (score_noDP == -1) {
                GET_TIME;
                int score_ODPD = ODPD(curr_read, curr_nbr, min, size_neighbour_in_symbols);
                STORE_ODPD_TIME;
                score = score_ODPD;
            }
            if (score <= min) {
                if (score < min) {
                    min = score;
                    nb_map = nb_map_start;
                }
                if (nb_map < MAX_DPU_RESULTS - 1) {
                    M_res[nb_map].num = M->num[current_read];
                    M_res[nb_map].coord = curr_coord;
                    M_res[nb_map].score = score;
                    nb_map++;
                    M_res[nb_map].num = -1;
                } else {
                    fprintf(stderr, "MAX_DPU_RESULTS reached!\n");
                    exit(-42);
                }
            }
        }
        current_read++;
    }
    WRITE_BACK_TIME(dpu_arg);
    return NULL;
}

void add_seed_to_simulation_requests(
    __attribute__((unused)) dispatch_request_t *requests, int num_read, int nb_read_written, index_seed_t *seed, int8_t *nbr)
{
    write_count(seed->num_dpu, nb_read_written, seed->nb_nbr);
    write_offset(seed->num_dpu, nb_read_written, seed->offset);
    write_num(seed->num_dpu, nb_read_written, num_read);
    write_num(seed->num_dpu, nb_read_written + 1, -1);
    write_neighbour_read(seed->num_dpu, nb_read_written, nbr);
}

void run_dpu_simulation(__attribute__((unused)) unsigned int dpu_offset, __attribute__((unused)) unsigned rank_id,
    __attribute__((unused)) unsigned int pass_id, int delta_neighbour, sem_t *dispatch_free_sem, sem_t *acc_wait_sem)
{
    unsigned int nb_dpu = get_nb_dpu();
    pthread_t thread_id[nb_dpu];
    align_on_dpu_arg_t thread_args[nb_dpu];

    sem_wait(acc_wait_sem);

    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        thread_args[numdpu].numdpu = numdpu;
        thread_args[numdpu].delta_neighbour = delta_neighbour;
#ifdef PRINT_NODP_ODPD_TIME
        thread_args[numdpu].nodp_time = 0.0;
        thread_args[numdpu].odpd_time = 0.0;
        thread_args[numdpu].nodp_call = 0;
        thread_args[numdpu].odpd_call = 0;
#endif
    }

    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        pthread_create(&thread_id[numdpu], NULL, align_on_dpu, &thread_args[numdpu]);
    }
    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        pthread_join(thread_id[numdpu], NULL);
    }

#ifdef PRINT_NODP_ODPD_TIME
    double odpd_total_time = 0.0;
    double nodp_total_time = 0.0;
    unsigned int nodp_total_call = 0;
    unsigned int odpd_total_call = 0;
    for (unsigned int numdpu = first_dpu; numdpu < nb_dpu; numdpu++) {
        odpd_total_time += thread_args[numdpu].odpd_time;
        nodp_total_time += thread_args[numdpu].nodp_time;
        nodp_total_call += thread_args[numdpu].nodp_call;
        odpd_total_call += thread_args[numdpu].odpd_call;
    }

    printf("nodp time: %lf sec (%u call)\n", nodp_total_time, nodp_total_call);
    printf("odpd time: %lf sec (%u call)\n", odpd_total_time, odpd_total_call);
#endif

    sem_post(dispatch_free_sem);
}

void init_backend_simulation() { malloc_dpu(get_nb_dpu()); }

void free_backend_simulation(unsigned int nb_dpu) { free_dpu(nb_dpu); }

void load_mram_simulation(__attribute__((unused)) unsigned int dpu_offset, __attribute__((unused)) unsigned int rank_id,
    __attribute__((unused)) int delta_neighbour)
{
    mem_dpu_t *MDPU;
    unsigned int nb_dpu = get_nb_dpu();

    for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
        MDPU = get_mem_dpu(each_dpu);
        mram_load((uint8_t *)(MDPU->neighbour_idx), each_dpu);
    }
}
