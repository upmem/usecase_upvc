#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

#include "upvc_dpu.h"
#include "code.h"
#include "index.h"
#include "genome.h"
#include "upvc.h"

#define MAX_SIZE_IDX_SEED (1000)

typedef struct seed_counter {
        int nb_seed;
        int seed_code;
} seed_counter_t;

static int cmp_seed_counter(void const *a, void const *b)
{
        seed_counter_t *seed_counter_a =  (seed_counter_t *)a;
        seed_counter_t *seed_counter_b =  (seed_counter_t *)b;
        if (seed_counter_a->nb_seed < seed_counter_b->nb_seed) {
                return 1;
        } else {
                return -1;
        }
}

index_seed_t **index_genome(genome_t *ref_genome, times_ctx_t *times_ctx, reads_info_t *reads_info)
{
        double t1,t2;
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        seed_counter_t *seed_counter;
        long dpu_workload[NB_DPU] = {0};
        int dpu_index_size[NB_DPU] = {0};
        int dpu_offset_in_memory[NB_DPU] = {0};
        int8_t buf_code_neighbour[size_neighbour];
        index_seed_t **index_seed = (index_seed_t **) malloc(sizeof(index_seed_t *)*NB_SEED);

        t1 = my_clock();

        /* Initialize the seed_counter table */
        seed_counter  = (seed_counter_t *)malloc(NB_SEED * sizeof(seed_counter_t));
        for (int i = 0; i < NB_SEED; i++) {
                seed_counter[i].nb_seed = 0;
                seed_counter[i].seed_code = i;
        }
        for (int i = 0; i < ref_genome->nb_seq; i++) {
                int sequence_start_idx = ref_genome->pt_seq[i];
                for (int sequence_idx = 0;
                     sequence_idx < ref_genome->len_seq[i] - size_neighbour - SIZE_SEED + 1;
                     sequence_idx++) {
                        int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                        if (seed_code >= 0) {
                                seed_counter[seed_code].nb_seed++;
                        }
                }
        }

        /* Create and initialize and link together all the seed */
        for (int i = 0; i < NB_SEED; i++) {
                int nb_index_needed = (seed_counter[i].nb_seed / MAX_SIZE_IDX_SEED) + 1;
                int nb_neighbour_per_index = (seed_counter[i].nb_seed / nb_index_needed) + 1;

                index_seed[i] = (index_seed_t *) malloc(sizeof(index_seed_t));
                index_seed[i]->nb_nbr = nb_neighbour_per_index;
                index_seed[i]->next = NULL;
                for (int j = 1; j < nb_index_needed; j++) {
                        index_seed_t *seed = (index_seed_t *) malloc(sizeof(index_seed_t));
                        seed->nb_nbr = nb_neighbour_per_index;
                        seed->next = index_seed[i];
                        index_seed[i] = seed;
                }
                index_seed[i]->nb_nbr = seed_counter[i].nb_seed - ((nb_index_needed - 1) * nb_neighbour_per_index);
        }

        /* Distribute indexs between DPUs */
        /* Sort the seeds from the most to the least used. */
        qsort(seed_counter, NB_SEED, sizeof(seed_counter_t), cmp_seed_counter);
        {
                int current_dpu = 0;
                for (int i = 0; i < NB_SEED; i++) {
                        int seed_code = seed_counter[i].seed_code;
                        int nb_seed_counted = seed_counter[i].nb_seed;
                        index_seed_t *seed = index_seed[seed_code];

                        while (seed != NULL) {
                                int dpu_for_current_seed = current_dpu;
                                if (nb_seed_counted != 0) {
                                        long max_dpu_workload = LONG_MAX;
                                        for (int j = 0; j < NB_DPU; j++) {
                                                if (dpu_workload[current_dpu] < max_dpu_workload) {
                                                        dpu_for_current_seed = current_dpu;
                                                        max_dpu_workload = dpu_workload[current_dpu];
                                                }
                                                current_dpu = (current_dpu + 1) % NB_DPU;
                                        }
                                }
                                dpu_index_size[dpu_for_current_seed] += seed->nb_nbr;
                                dpu_workload[dpu_for_current_seed] += (long) (seed->nb_nbr * nb_seed_counted);
                                seed->num_dpu = dpu_for_current_seed;
                                seed = seed->next;
                                current_dpu = (current_dpu + 1) % NB_DPU;
                        }
                }
        }

        /* Compute offset in the DPUs memories */
        for (int i = 0; i < NB_SEED; i++) {
                index_seed_t *seed = index_seed[i];
                while (seed != NULL) {
                        seed->offset = dpu_offset_in_memory[seed->num_dpu];
                        dpu_offset_in_memory[seed->num_dpu] += seed->nb_nbr;
                        seed = seed->next;
                }
        }

        malloc_dpu(reads_info);
        for (int i = 0; i < NB_DPU; i++) {
                malloc_neighbour_idx(i, dpu_index_size[i], reads_info);
        }

        /* Writing data in DPUs memories */
        memset(seed_counter, 0, sizeof(seed_counter_t) * NB_SEED);
        for (int seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
                int sequence_start_idx = ref_genome->pt_seq[seq_number];
                for (int sequence_idx = 0;
                     sequence_idx < ref_genome->len_seq[seq_number] - size_neighbour - SIZE_SEED + 1;
                     sequence_idx++) {
                        index_seed_t *seed;
                        int total_nb_neighbour = 0;
                        int align_idx;
                        int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);

                        if (seed_code < 0) {
                                continue;
                        }

                        seed = index_seed[seed_code];
                        while (seed != NULL) {
                                if (seed_counter[seed_code].nb_seed < seed->nb_nbr + total_nb_neighbour)
                                        break;
                                total_nb_neighbour += seed->nb_nbr;
                                seed = seed->next;
                        }
                        align_idx = seed->offset + seed_counter[seed_code].nb_seed - total_nb_neighbour;

                        /* TODO: to be replace by transfer with DPUs memory */
                        code_neighbour(&ref_genome->data[sequence_start_idx+sequence_idx+SIZE_SEED],
                                       buf_code_neighbour,
                                       reads_info);
                        write_neighbour_idx(seed->num_dpu, align_idx, buf_code_neighbour, reads_info);
                        write_coordinate(seed->num_dpu, align_idx, ( ((long)seq_number) << 32) + sequence_idx);
                        seed_counter[seed_code].nb_seed++;
                }
        }

        free(seed_counter);

        t2 = my_clock();
        times_ctx->index_genome = t2 - t1;

        return index_seed;
}

void free_index(index_seed_t **index_seed)
{
        for (int i = 0; i < NB_SEED; i++) {
                index_seed_t *seed = index_seed[i];
                while (seed != NULL) {
                        index_seed_t *next_seed = seed->next;
                        free(seed);
                        seed = next_seed;
                }
        }
        free(index_seed);
}
