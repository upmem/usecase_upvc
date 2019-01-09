/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "upvc_dpu.h"
#include "code.h"
#include "index.h"
#include "genome.h"
#include "upvc.h"
#include "vmi.h"
#include "mram_dpu.h"

#include "common.h"

#define MAX_SIZE_IDX_SEED (1000)
#define SEED_FILE ("seeds.txt")

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

index_seed_t **load_index_seeds()
{
        FILE *f = fopen(SEED_FILE, "r");
        index_seed_t **index_seed;
        unsigned int max_dpu = 0;
        double t1, t2;

        assert(f != NULL);

        printf("Loading index seeds\n");
        t1 = my_clock();

        index_seed = (index_seed_t **) calloc(NB_SEED, sizeof(index_seed_t *));

        { /* First line is just a comment, skip */
                char line[512];
                assert(fgets(line, sizeof(line), f) != NULL);
        }

        {
                unsigned int seed_id = 0, dpu = 0, offset = 0, nb_nbr = 0;
                while (fscanf(f, "%u %u %u %u", &seed_id, &dpu, &offset, &nb_nbr) == 4) {
                        index_seed_t *seed = index_seed[seed_id];
                        if (seed == NULL) {
                                index_seed[seed_id] = seed = (index_seed_t *) malloc(sizeof(index_seed_t));
                        } else {
                                while (seed->next != NULL) {
                                        seed = seed->next;
                                }
                                seed->next = (index_seed_t *) malloc(sizeof(index_seed_t));
                                seed = seed->next;
                        }

                        if (dpu > max_dpu) {
                                max_dpu = dpu;
                        }
                        seed->num_dpu = dpu;
                        seed->nb_nbr = nb_nbr;
                        seed->offset = offset;
                        seed->next = NULL;
                }
        }

        set_nb_dpu(max_dpu + 1);
        if (get_simulation_mode()) {
                set_nb_dpus_per_run(max_dpu + 1);
        }
        printf(" - nb_dpu: %u\n", max_dpu + 1);

        fclose(f);

        t2 = my_clock();
        printf(" - time: %lf\n", t2 - t1);

        return index_seed;
}

void save_index_seeds(index_seed_t **index_seed)
{
    index_seed_t *seed;
    FILE *f = fopen(SEED_FILE, "w");
    assert(f != NULL);

    fprintf(f, "# seed dpu offset nb_nbr\n");
    for (int i = 0; i < NB_SEED; i++) {
        seed = index_seed[i];
        while (seed != NULL) {
            fprintf(f, "%u %u %u %u\n", i, seed->num_dpu, seed->offset, seed->nb_nbr);
            seed = seed->next;
        }
    }

    fclose(f);
}

static vmi_t *init_vmis(unsigned int nb_dpu)
{
        vmi_t *vmis = (vmi_t *) calloc(nb_dpu, sizeof(vmi_t));
        char vmi_name[8];
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                (void) sprintf(vmi_name, "%04u", dpuno);
                vmi_create(vmi_name, vmis + dpuno);
        }
        return vmis;
}

static void dump_mdpu_images_into_mram_files(vmi_t *vmis, unsigned int *nbr_counts, unsigned int nb_dpu, reads_info_t *reads_info)
{
        mram_info_t *mram_image = (mram_info_t *)malloc(MRAM_SIZE);
        printf("Creating MRAM images\n");
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
                vmi_t *this_vmi = vmis + each_dpu;
                mram_copy_vmi(mram_image, this_vmi, nbr_counts[each_dpu], reads_info);
                mram_save(mram_image, each_dpu);
        }
        free(mram_image);
}

static void free_vmis(vmi_t *vmis, unsigned int nb_dpu, unsigned int *nb_neighbours, reads_info_t *reads_info)
{
        dump_mdpu_images_into_mram_files(vmis, nb_neighbours, nb_dpu, reads_info);
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                vmi_delete(vmis + dpuno);
        }
        free(vmis);
}

static void write_vmi(vmi_t *vmis,
                          unsigned int dpuno,
                          unsigned int k,
                          int8_t *nbr,
                          dpu_result_coord_t coord,
                          reads_info_t *reads_info)
{
        unsigned int size_neighbour_in_bytes = reads_info->size_neighbour_in_bytes;
        unsigned int out_len = ALIGN_DPU(sizeof(dpu_result_coord_t) + size_neighbour_in_bytes);
        uint64_t temp_buff[out_len / sizeof(uint64_t)];
        memset(temp_buff, 0, out_len);
        temp_buff[0] = coord.coord;
        memcpy(&temp_buff[1], nbr, (size_t) size_neighbour_in_bytes);
        vmi_write(vmis + dpuno, k * out_len, temp_buff, out_len);
}

index_seed_t **index_genome(genome_t *ref_genome,
                            int nb_dpu,
                            times_ctx_t *times_ctx,
                            reads_info_t *reads_info)
{
        double t1,t2;
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        seed_counter_t *seed_counter;
        long dpu_workload[nb_dpu];
        int dpu_index_size[nb_dpu];
        int dpu_offset_in_memory[nb_dpu];
        int8_t buf_code_neighbour[size_neighbour];
        index_seed_t **index_seed = (index_seed_t **) malloc(sizeof(index_seed_t *)*NB_SEED);
        vmi_t *vmis = NULL;

        printf("Index genome\n");
        t1 = my_clock();

        vmis = init_vmis(nb_dpu);

        memset(&dpu_workload, 0, nb_dpu * sizeof(long));
        memset(&dpu_index_size, 0, nb_dpu * sizeof(int));
        memset(&dpu_offset_in_memory, 0, nb_dpu * sizeof(int));

        /* Initialize the seed_counter table */
        printf("\tInitialize the seed_counter table\n");
        seed_counter  = (seed_counter_t *)malloc(NB_SEED * sizeof(seed_counter_t));
        for (int i = 0; i < NB_SEED; i++) {
                seed_counter[i].nb_seed = 0;
                seed_counter[i].seed_code = i;
        }
        for (uint32_t i = 0; i < ref_genome->nb_seq; i++) {
                uint64_t sequence_start_idx = ref_genome->pt_seq[i];
                for (uint64_t sequence_idx = 0;
                     sequence_idx < ref_genome->len_seq[i] - size_neighbour - SIZE_SEED + 1;
                     sequence_idx++) {
                        int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                        if (seed_code >= 0) {
                                seed_counter[seed_code].nb_seed++;
                        }
                }
        }

        /* Create and initialize and link together all the seed */
        printf("\tCreate, initialize and link together all the seed\n");
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
        printf("\tDistribute indexs between DPUs\n");
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
                                        for (int j = 0; j < nb_dpu; j++) {
                                                if (dpu_workload[current_dpu] < max_dpu_workload) {
                                                        dpu_for_current_seed = current_dpu;
                                                        max_dpu_workload = dpu_workload[current_dpu];
                                                }
                                                current_dpu = (current_dpu + 1) % nb_dpu;
                                        }
                                }
                                dpu_index_size[dpu_for_current_seed] += seed->nb_nbr;
                                dpu_workload[dpu_for_current_seed] += (long) (seed->nb_nbr * nb_seed_counted);
                                seed->num_dpu = dpu_for_current_seed;
                                seed = seed->next;
                                current_dpu = (current_dpu + 1) % nb_dpu;
                        }
                }
        }

        /* Compute offset in the DPUs memories */
        printf("\tCompute offset in the DPUs memories\n");
        for (int i = 0; i < NB_SEED; i++) {
                index_seed_t *seed = index_seed[i];
                while (seed != NULL) {
                        seed->offset = dpu_offset_in_memory[seed->num_dpu];
                        dpu_offset_in_memory[seed->num_dpu] += seed->nb_nbr;
                        seed = seed->next;
                }
        }

        unsigned int nb_neighbours[nb_dpu];
        memset(nb_neighbours, 0, nb_dpu * sizeof(unsigned int));

        /* Writing data in DPUs memories */
        printf("\tWriting data in DPUs memories\n");
        memset(seed_counter, 0, sizeof(seed_counter_t) * NB_SEED);
        for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
                uint64_t sequence_start_idx = ref_genome->pt_seq[seq_number];
                for (uint64_t sequence_idx = 0;
                     sequence_idx < ref_genome->len_seq[seq_number] - size_neighbour - SIZE_SEED + 1;
                     sequence_idx++) {
                        index_seed_t *seed;
                        int total_nb_neighbour = 0;
                        int align_idx;
                        int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                        dpu_result_coord_t coord_var =
                                {
                                 .seq_nr = seq_number,
                                 .seed_nr = sequence_idx,
                                };

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

                        code_neighbour(&ref_genome->data[sequence_start_idx+sequence_idx+SIZE_SEED],
                                       buf_code_neighbour,
                                       reads_info);
                        if (sequence_idx % 1000 == 0)
                                printf("\r\t%lli/%lli %i/%i         ",
                                       (unsigned long long)sequence_idx,
                                       (unsigned long long)ref_genome->len_seq[seq_number] - size_neighbour - SIZE_SEED + 1,
                                       seq_number,
                                       ref_genome->nb_seq);
                        write_vmi(vmis,
                                  seed->num_dpu,
                                  align_idx,
                                  buf_code_neighbour,
                                  coord_var,
                                  reads_info);
                        nb_neighbours[seed->num_dpu]++;

                        seed_counter[seed_code].nb_seed++;
                }
                printf("\n");
        }

        free_vmis(vmis, nb_dpu, nb_neighbours, reads_info);
        free(seed_counter);

        t2 = my_clock();
        times_ctx->index_genome = t2 - t1;

        printf(" - time: %lf\n", times_ctx->index_genome);

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

void print_index_seeds(index_seed_t **index_seed, FILE *out, reads_info_t *reads_info)
{
    for (int i = 0; i < NB_SEED; i++) {
        fprintf(out, "SEED=%u\n", i);
        index_seed_t *seed = index_seed[i];
        while (seed != NULL) {
            fprintf(out, "\tPU %u @%u[%u]\n", seed->num_dpu, seed->offset, seed->nb_nbr);
            print_neighbour_idx(seed->num_dpu, seed->offset, seed->nb_nbr, out, reads_info);
            fprintf(out, "\t");
            print_coordinates(seed->num_dpu, seed->offset, seed->nb_nbr, out, reads_info);
            seed = seed->next;
        }
    }
}
