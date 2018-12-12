#include <stdint.h>
#include <stdio.h>
#include <pthread.h>

#include "dpus_mgmt.h"
#include "upvc_dpu.h"
#include "genome.h"
#include "upvc.h"
#include "vmi.h"
#include "dispatch.h"
#include "index.h"
#include "parse_args.h"
#include "backends_functions.h"

#include "common.h"

#define MAX_SCORE        40

/**
 * @brief Arguments needed for the thread that will be offload onto a DPU.
 */
typedef struct align_on_dpu_arg {
        reads_info_t reads_info;
        int numdpu;
} align_on_dpu_arg_t;

static int min (int a, int b)
{
        return a < b ? a : b;
}

#define PQD_INIT_VAL (99)
static void ODPD_compute(int i,
                         int j,
                         int8_t *s1,
                         int8_t *s2,
                         int *Pppj,
                         int Pppjm,
                         int *Qppj,
                         int Qlpj,
                         int *Dppj,
                         int Dppjm,
                         int Dlpj,
                         int Dlpjm)
{
        int d = Dlpjm;
        int QP;
        *Pppj = min(Dppjm + COST_GAPO, Pppjm + COST_GAPE);
        *Qppj = min(Dlpj  + COST_GAPO, Qlpj  + COST_GAPE);
        QP = min(*Pppj, *Qppj);
        if ( ( (s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3)
             !=
             ( (s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3) ) {
                d += COST_SUB;
        }
        *Dppj = min(d, QP);
}

/**
 * @brief Compute the alignment distance by dynamical programming on the diagonals of the matrix.
 * Stops when score is greater than max_score
 */
static int ODPD(int8_t *s1, int8_t *s2, int max_score, reads_info_t *reads_info)
{
        int size_neighbour = reads_info->size_neighbour_in_32bits_words;
        int matrix_size = size_neighbour+1;
        int D[2][matrix_size];
        int P[2][matrix_size];
        int Q[2][matrix_size];
        int diagonal = (NB_DIAG / 2) + 1;

        for (int j = 0; j <= diagonal; j++) {
                P[0][j] = PQD_INIT_VAL;
                Q[0][j] = PQD_INIT_VAL;
                D[0][j] = j * COST_SUB;
        }
        P[1][0]=PQD_INIT_VAL;
        Q[1][0]=PQD_INIT_VAL;

        for (int i = 1; i < diagonal; i++) {
                int min_score = PQD_INIT_VAL;
                int pp = i % 2;
                int lp = (i - 1) % 2;
                D[pp][0] = i * COST_SUB;
                for (int j = 1; j < i + diagonal; j++) {
                        ODPD_compute(i, j, s1, s2,
                                     &P[pp][j], P[pp][j - 1],
                                     &Q[pp][j], Q[lp][j],
                                     &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
                        if (D[pp][j] < min_score) {
                                min_score = D[pp][j];
                        }
                }
                Q[pp][i + diagonal]=PQD_INIT_VAL;
                D[pp][i + diagonal]=PQD_INIT_VAL;
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
                        ODPD_compute(i, j, s1, s2,
                                     &P[pp][j], P[pp][j - 1],
                                     &Q[pp][j], Q[lp][j],
                                     &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
                        if (D[pp][j] < min_score) {
                                min_score = D[pp][j];
                        }
                }
                Q[pp][i + diagonal]=PQD_INIT_VAL;
                D[pp][i + diagonal]=PQD_INIT_VAL;
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
                        ODPD_compute(i, j, s1, s2,
                                     &P[pp][j], P[pp][j - 1],
                                     &Q[pp][j], Q[lp][j],
                                     &D[pp][j], D[pp][j - 1], D[lp][j], D[lp][j - 1]);
                }
                if (D[pp][matrix_size - 1] < min_score)
                        min_score = D[pp][matrix_size - 1];
        }
        int i=matrix_size - 1;
        int pp = i%2;
        for (int j = i + 1 - diagonal; j < matrix_size; j++) {
                if (D[pp][j] < min_score) {
                        min_score = D[pp][j];
                }
        }

        return min_score;
}

static int translation_table[256] = {0 , 10 , 10 , 10 , 10 , 20 , 20 , 20 , 10 , 20 , 20 , 20 , 10 , 20 , 20 , 20 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
                                     20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 };

/**
 * @brief Optimized version of ODPD (if no INDELS)
 * If it detects INDELS, return -1. In this case we will need the run the full ODPD.
 */
static int noDP(int8_t *s1, int8_t *s2, int max_score, reads_info_t *reads_info)
{
        int score = 0;
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        int delta_neighbour = reads_info->delta_neighbour_in_bytes;
        for (int i = 0; i < size_neighbour - delta_neighbour; i++) {
                int s_xor = ((int) (s1[i] ^ s2[i])) & 0xFF;
                int s_translated = translation_table[s_xor];
                if (s_translated > COST_SUB) {
                        int j = i + 1;
                        /* INDELS detection */
                        if (j < size_neighbour - delta_neighbour - 3) {
                                int s1_val = ((int *) (&s1[j]))[0];
                                int s2_val = ((int *) (&s2[j]))[0];
                                if ( ((s1_val ^ (s2_val >> 2)) & 0x3FFFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s1_val ^ (s2_val >> 4)) & 0xFFFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s1_val ^ (s2_val >> 6)) & 0x3FFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s1_val ^ (s2_val >> 8)) & 0xFFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s2_val ^ (s1_val >> 2)) & 0x3FFFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s2_val ^ (s1_val >> 4)) & 0xFFFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s2_val ^ (s1_val >> 6)) & 0x3FFFFFF) == 0) {
                                        return -1;
                                }
                                if ( ((s2_val ^ (s1_val >> 8)) & 0xFFFFFF) == 0) {
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
        reads_info_t *reads_info = &(dpu_arg->reads_info);
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        mem_dpu_t *M = get_mem_dpu(numdpu);
        dpu_result_out_t *M_res = get_mem_dpu_res(numdpu);
        int current_read = 0;
        /* int debug_nb = -1; */

        M_res[nb_map].num = -1;
        while (M->num[current_read] != -1) {
                int offset = M->offset[current_read]; /* address of the first neighbour */
                int min = MAX_SCORE;
                int nb_map_start = nb_map;
                /* print_neighbour((uint8_t*)&M->neighbour_read[current_read*size_neighbour], stdout, reads_info); */
                /* debug_nb++; */
                /* if (toto++ > 5) */
                /*         exit(0); */
                for (int nb_neighbour = 0; nb_neighbour < M->count[current_read]; nb_neighbour++) {
                        /* print_neighbour((uint8_t*)&M->neighbour_idx[(offset+nb_neighbour)*size_neighbour], stdout, reads_info); */
                        int score_noDP =noDP(&M->neighbour_read[current_read*size_neighbour],
                                             &M->neighbour_idx[(offset+nb_neighbour)*size_neighbour],
                                             min,
                                             reads_info);
                        /* printf("noDP score = %i\n", score_noDP); */
                        int score = score_noDP;
                        if (score_noDP == -1) {
                                int score_ODPD = ODPD(&M->neighbour_read[current_read*size_neighbour],
                                                      &M->neighbour_idx[(offset+nb_neighbour)*size_neighbour],
                                                      min,
                                                      reads_info);
                                score = score_ODPD;
                                /* printf("ODPD score = %i\n", score_ODPD); */
                        }
                        if (score <= min) {
                                if (score <  min) {
                                        /* printf("mini\n"); */
                                        min = score;
                                        nb_map = nb_map_start;
                                }
                                if (nb_map < MAX_DPU_RESULTS-1) {
                                        /* printf("[0] nr=%u ix=%u num=%u ", debug_nb, nb_neighbour, M->num[current_read]); */
                                        /* printf("offset=%u seed=%li seq=%li score=%i\n", */
                                        /*        offset, */
                                        /*        M->coordinate[offset+nb_neighbour] &0xffffffff, */
                                        /*        (M->coordinate[offset+nb_neighbour]>>32) &0xffffffff, */
                                        /*        score); */
                                        M_res[nb_map].num = M->num[current_read];
                                        M_res[nb_map].coord.coord = M->coordinate[offset+nb_neighbour];
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
        return NULL;
}


vmi_t *init_vmis_simulation(__attribute__((unused)) unsigned int nb_dpu) { return NULL; }

void free_vmis_simulation(__attribute__((unused)) vmi_t *vmis,
                          __attribute__((unused)) unsigned int nb_dpu,
                          __attribute__((unused)) unsigned int *nb_neighbours,
                          __attribute__((unused)) reads_info_t *reads_info) {}

void write_vmi_simulation(__attribute__((unused)) vmi_t *vmis,
                          unsigned int dpuno,
                          unsigned int align_idx,
                          int8_t *nbr,
                          uint64_t coords,
                          reads_info_t *reads_info)
{
        write_neighbour_idx(dpuno, align_idx, nbr, reads_info);
        write_coordinate(dpuno, align_idx, coords);
}

index_seed_t **get_index_seed_simulation(unsigned int nb_dpu,
                                         genome_t *ref_genome,
                                         reads_info_t *reads_info,
                                         times_ctx_t *times_ctx,
                                         backends_functions_t *backends_functions)
{
        return index_genome(ref_genome, nb_dpu, times_ctx, reads_info, backends_functions);
}

void add_seed_to_simulation_requests(__attribute__((unused)) dispatch_request_t *requests,
                                     __attribute__((unused)) int num_dpu,
                                     int num_read,
                                     int nb_read_written,
                                     index_seed_t *seed,
                                     int8_t *nbr,
                                     reads_info_t *reads_info)
{
        write_count(seed->num_dpu, nb_read_written, seed->nb_nbr);
        write_offset(seed->num_dpu, nb_read_written, seed->offset);
        write_num(seed->num_dpu, nb_read_written, num_read);
        write_neighbour_read(seed->num_dpu, nb_read_written, nbr, reads_info);
}

void run_dpu_simulation(__attribute__((unused)) dispatch_t dispatch,
                        __attribute__((unused)) devices_t *devices,
                        unsigned int nb_dpu,
                        times_ctx_t *times_ctx,
                        reads_info_t *reads_info)
{
        double t1, t2;
        pthread_t thread_id[nb_dpu];
        align_on_dpu_arg_t thread_args[nb_dpu];
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        unsigned int first_dpu = ((DEBUG_FIRST_RUN != -1) ? (DEBUG_FIRST_RUN * nb_dpus_per_run) : 0);

        t1 = my_clock();

        for (unsigned int numdpu = first_dpu; numdpu < nb_dpu; numdpu++) {
                thread_args[numdpu].reads_info = *reads_info;
                thread_args[numdpu].numdpu = numdpu;
        }

        for (unsigned int numdpu = first_dpu; numdpu < nb_dpu; numdpu++) {
                pthread_create(&thread_id[numdpu], NULL, align_on_dpu, &thread_args[numdpu]);
        }
        for (unsigned int numdpu = first_dpu; numdpu < nb_dpu; numdpu++) {
                pthread_join(thread_id[numdpu], NULL);
        }

        for (; first_dpu < nb_dpu; first_dpu += nb_dpus_per_run) {
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        printf("LOG DPU=%u\n", this_dpu);
                        int k = 0;
                        int num_read = read_out_num(this_dpu, k);
                        while (num_read != -1) {
                                printf("R: %u %u %lu\n", num_read, read_out_score(this_dpu, k), read_out_coord(this_dpu, k).coord);
                                k++;
                                num_read = read_out_num(this_dpu, k);
                        }
                }
        }

        t2 = my_clock();
        times_ctx->map_read = t2 - t1;
        times_ctx->tot_map_read += t2 - t1;
}

devices_t *init_devices_simulation(__attribute__((unused)) unsigned int nb_dpu,
                                   __attribute__((unused)) const char *dpu_binary)
{
        return NULL;
}

void free_devices_simulation(__attribute__((unused)) devices_t *devices) {}
