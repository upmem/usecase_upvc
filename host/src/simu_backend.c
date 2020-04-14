/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>
#include <stdio.h>

#include "accumulateread.h"
#include "dispatch.h"
#include "index.h"
#include "mram_dpu.h"
#include "simu_backend.h"
#include "upvc.h"

#include <dpu.h>

#include "common.h"

#define MAX_SCORE 40

static coords_and_nbr_t **mrams;
static int delta_neighbour;
static int pass_id;

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
static int noDP(int8_t *s1, int8_t *s2, int max_score)
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
    int numdpu = *(int *)arg;
    int size_neighbour_in_symbols = SIZE_IN_SYMBOLS(delta_neighbour);
    dispatch_request_t *requests = dispatch_get(numdpu, pass_id);
    acc_results_t *acc_res = accumulate_get_buffer(numdpu, pass_id);

    for (unsigned int each_request_read = 0; each_request_read < requests->nb_reads; each_request_read++) {
        dpu_request_t *curr_request = &(requests->dpu_requests[each_request_read]);
        int min = MAX_SCORE;
        int nb_map_start = nb_map;
        int8_t *curr_read = (int8_t *)&curr_request->nbr[0];
        for (unsigned int nb_neighbour = 0; nb_neighbour < curr_request->count; nb_neighbour++) {
            coords_and_nbr_t *coord_and_nbr = &(mrams[numdpu][curr_request->offset + nb_neighbour]);
            int8_t *curr_nbr = (int8_t *)&coord_and_nbr->nbr[0];

            int score = noDP(curr_read, curr_nbr, min);
            if (score == -1) {
                score = ODPD(curr_read, curr_nbr, min, size_neighbour_in_symbols);
            }
            if (score > min)
                continue;

            if (score < min) {
                min = score;
                nb_map = nb_map_start;
            }

            if (nb_map >= MAX_DPU_RESULTS - 1) {
                ERROR_EXIT(-42, "%s:[P%u, DPU#%u]: MAX_DPU_RESULTS reached!", __func__, pass_id, numdpu);
            }

            dpu_result_out_t *result = &acc_res->results[nb_map++];
            result->num = curr_request->num;
            result->coord = coord_and_nbr->coord;
            result->score = score;
        }
    }

    acc_res->results[nb_map].num = -1;
    acc_res->nb_res = nb_map;

    return NULL;
}

void run_dpu_simulation(__attribute__((unused)) unsigned int dpu_offset, __attribute__((unused)) unsigned rank_id,
    unsigned int _pass_id, sem_t *dispatch_free_sem, sem_t *acc_wait_sem)
{
    unsigned int nb_dpu = index_get_nb_dpu();
    pthread_t thread_id[nb_dpu];
    int thread_args[nb_dpu];
    pass_id = _pass_id;

    sem_wait(acc_wait_sem);

    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        thread_args[numdpu] = numdpu;
    }

    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        pthread_create(&thread_id[numdpu], NULL, align_on_dpu, &thread_args[numdpu]);
    }
    for (unsigned int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        pthread_join(thread_id[numdpu], NULL);
    }

    sem_post(dispatch_free_sem);
}

void init_backend_simulation(unsigned int *nb_dpus_per_run, unsigned int *nb_ranks_per_run)
{
    *nb_ranks_per_run = 1;
    *nb_dpus_per_run = index_get_nb_dpu();
    mrams = (coords_and_nbr_t **)malloc(index_get_nb_dpu() * sizeof(coords_and_nbr_t *));
    assert(mrams != NULL);
}

void free_backend_simulation()
{
    for (unsigned int each_dpu = 0; each_dpu < index_get_nb_dpu(); each_dpu++) {
        free(mrams[each_dpu]);
    }
    free(mrams);
}

void load_mram_simulation(
    __attribute__((unused)) unsigned int dpu_offset, __attribute__((unused)) unsigned int rank_id, int _delta_neighbour)
{
    delta_neighbour = _delta_neighbour;
    for (unsigned int each_dpu = 0; each_dpu < index_get_nb_dpu(); each_dpu++) {
        mram_load((uint8_t **)&mrams[each_dpu], each_dpu);
    }
}
