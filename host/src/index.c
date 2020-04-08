/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "index.h"
#include "code.h"
#include "genome.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "upvc.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"

#define NB_SEED (1 << (SIZE_SEED << 1)) /* NB_SEED = 4 ^ (SIZE_SEED) */

#define MAX_SIZE_IDX_SEED (1000)
#define SEED_FILE ("seeds.txt")

typedef int vmi_t;

static index_seed_t **index_seed;
static vmi_t *vmis = NULL;

index_seed_t **index_get() { return index_seed; }

typedef struct seed_counter {
    int nb_seed;
    int seed_code;
} seed_counter_t;

static int cmp_seed_counter(void const *a, void const *b)
{
    seed_counter_t *seed_counter_a = (seed_counter_t *)a;
    seed_counter_t *seed_counter_b = (seed_counter_t *)b;
    if (seed_counter_a->nb_seed < seed_counter_b->nb_seed) {
        return 1;
    } else {
        return -1;
    }
}

void index_load()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);
    FILE *f = fopen(SEED_FILE, "r");
    assert(f != NULL);

    unsigned int max_dpu = 0;
    uint32_t seed_id = 0;

    index_seed = (index_seed_t **)calloc(NB_SEED, sizeof(index_seed_t *));
    assert(index_seed != NULL);

    index_seed_t *new_seed = (index_seed_t *)malloc(sizeof(index_seed_t));
    assert(new_seed != NULL);
    while (fread(&seed_id, sizeof(uint32_t), 1, f) == 1) {
        size_t read_nmemb = fread(new_seed, sizeof(uint32_t), 3, f);
        assert(read_nmemb == 3);
        new_seed->next = index_seed[seed_id];
        index_seed[seed_id] = new_seed;

        if (new_seed->num_dpu > max_dpu) {
            max_dpu = new_seed->num_dpu;
        }
        new_seed = (index_seed_t *)malloc(sizeof(index_seed_t));
        assert(new_seed != NULL);
    }
    free(new_seed);

    set_nb_dpu(max_dpu + 1);
    printf("\tnb_dpu: %u\n", max_dpu + 1);

    fclose(f);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void index_save()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);
    FILE *f = fopen(SEED_FILE, "w");
    assert(f != NULL);

    for (unsigned int i = 0; i < NB_SEED; i++) {
        index_seed_t *seed = index_seed[i];
        while (seed != NULL) {
            fwrite(&i, sizeof(uint32_t), 1, f);
            fwrite(seed, sizeof(uint32_t), 3, f);
            seed = seed->next;
        }
    }

    fclose(f);
    printf("\ttime: %lf\n", my_clock() - start_time);
}

static void init_vmis(unsigned int nb_dpu)
{
    vmis = (vmi_t *)calloc(nb_dpu, sizeof(vmi_t));
    assert(vmis != NULL);
    for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
        char file_name[FILE_NAME_SIZE];
        vmis[dpuno] = open(make_mram_file_name(file_name, dpuno), O_WRONLY | O_CREAT, S_IRUSR | S_IRGRP | S_IROTH);
        if (vmis[dpuno] == -1) {
            ERROR_EXIT(-1, "%s: open failed (%s)", file_name, strerror(errno));
        }
    }
}

static void free_vmis(unsigned int nb_dpu)
{
    for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
        close(vmis[dpuno]);
    }
    free(vmis);
}

static void write_vmi(unsigned int dpuno, unsigned int k, int8_t *nbr, dpu_result_coord_t coord)
{
    unsigned int out_len = ALIGN_DPU(sizeof(dpu_result_coord_t) + SIZE_NEIGHBOUR_IN_BYTES);
    uint64_t temp_buff[out_len / sizeof(uint64_t)];
    memset(temp_buff, 0, out_len);
    temp_buff[0] = coord.coord;
    memcpy(&temp_buff[1], nbr, SIZE_NEIGHBOUR_IN_BYTES);

    vmi_t vmi = vmis[dpuno];
    lseek(vmi, k * out_len, SEEK_SET);
    ssize_t write_size = write(vmi, temp_buff, out_len);
    if (write_size != out_len) {
        char file_name[FILE_NAME_SIZE];
        ERROR_EXIT(-1, "%s: did not write the expected number of bytes (expected %u got %li, errno: %s)",
            make_mram_file_name(file_name, dpuno), out_len, write_size, strerror(errno));
    }
}

void index_init(int nb_dpu)
{
    double start_time = my_clock();
    printf("%s(%i):\n", __func__, nb_dpu);
    seed_counter_t *seed_counter;
    int8_t buf_code_neighbour[SIZE_NEIGHBOUR_IN_BYTES];

    index_seed = (index_seed_t **)malloc(sizeof(index_seed_t *) * NB_SEED);
    assert(index_seed != NULL);

    init_vmis(nb_dpu);

    /* Initialize the seed_counter table */
    {
        double init_seed_counter_time = my_clock();
        printf("\tInitialize the seed_counter table\n");
        seed_counter = (seed_counter_t *)malloc(NB_SEED * sizeof(seed_counter_t));
        for (int i = 0; i < NB_SEED; i++) {
            seed_counter[i].nb_seed = 0;
            seed_counter[i].seed_code = i;
        }

        genome_t *ref_genome = genome_get();
        for (uint32_t i = 0; i < ref_genome->nb_seq; i++) {
            uint64_t sequence_start_idx = ref_genome->pt_seq[i];
            for (uint64_t sequence_idx = 0; sequence_idx < ref_genome->len_seq[i] - SIZE_NEIGHBOUR_IN_BYTES - SIZE_SEED + 1;
                 sequence_idx++) {
                int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                if (seed_code >= 0) {
                    seed_counter[seed_code].nb_seed++;
                }
            }
        }
        printf("\t\ttime: %lf s\n", my_clock() - init_seed_counter_time);
    }

    /* Create and initialize and link together all the seed */
    {
        double create_init_link_all_seed_time = my_clock();
        printf("\tCreate, initialize and link together all the seed\n");

        for (int i = 0; i < NB_SEED; i++) {
            int nb_index_needed = (seed_counter[i].nb_seed / MAX_SIZE_IDX_SEED) + 1;
            int nb_neighbour_per_index = (seed_counter[i].nb_seed / nb_index_needed) + 1;

            index_seed[i] = (index_seed_t *)malloc(sizeof(index_seed_t));
            assert(index_seed[i] != NULL);
            index_seed[i]->nb_nbr = nb_neighbour_per_index;
            index_seed[i]->next = NULL;
            for (int j = 1; j < nb_index_needed; j++) {
                index_seed_t *seed = (index_seed_t *)malloc(sizeof(index_seed_t));
                assert(seed != NULL);
                seed->nb_nbr = nb_neighbour_per_index;
                seed->next = index_seed[i];
                index_seed[i] = seed;
            }
            index_seed[i]->nb_nbr = seed_counter[i].nb_seed - ((nb_index_needed - 1) * nb_neighbour_per_index);
        }
        printf("\t\ttime: %lf s\n", my_clock() - create_init_link_all_seed_time);
    }

    /* Distribute indexs between DPUs */
    /* Sort the seeds from the most to the least used. */
    {
        double distribute_index_time = my_clock();
        printf("\tDistribute indexs between DPUs\n");

        qsort(seed_counter, NB_SEED, sizeof(seed_counter_t), cmp_seed_counter);

        long *dpu_workload = (long *)calloc(nb_dpu, sizeof(long));
        int *dpu_index_size = (int *)calloc(nb_dpu, sizeof(int));
        assert(dpu_workload != NULL && dpu_index_size != NULL);

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
                dpu_workload[dpu_for_current_seed] += (long)(seed->nb_nbr * nb_seed_counted);
                seed->num_dpu = dpu_for_current_seed;
                seed = seed->next;
                current_dpu = (current_dpu + 1) % nb_dpu;
            }
        }
        free(dpu_workload);
        free(dpu_index_size);
        printf("\t\ttime: %lf s\n", my_clock() - distribute_index_time);
    }

    /* Compute offset in the DPUs memories */
    {
        double offset_dpu_memories_time = my_clock();
        printf("\tCompute offset in the DPUs memories\n");
        int *dpu_offset_in_memory;
        dpu_offset_in_memory = (int *)calloc(nb_dpu, sizeof(int));
        for (int i = 0; i < NB_SEED; i++) {
            index_seed_t *seed = index_seed[i];
            while (seed != NULL) {
                seed->offset = dpu_offset_in_memory[seed->num_dpu];
                dpu_offset_in_memory[seed->num_dpu] += seed->nb_nbr;
                seed = seed->next;
            }
        }
        free(dpu_offset_in_memory);
        printf("\t\ttime: %lf s\n", my_clock() - offset_dpu_memories_time);
    }

    /* Writing data in DPUs memories */
    {
        double write_in_memories_time = my_clock();
        printf("\tWriting data in DPUs memories\n");

        memset(seed_counter, 0, sizeof(seed_counter_t) * NB_SEED);
        genome_t *ref_genome = genome_get();
        for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
            uint64_t sequence_start_idx = ref_genome->pt_seq[seq_number];
            for (uint64_t sequence_idx = 0;
                 sequence_idx < ref_genome->len_seq[seq_number] - SIZE_NEIGHBOUR_IN_BYTES - SIZE_SEED + 1; sequence_idx++) {
                index_seed_t *seed;
                int total_nb_neighbour = 0;
                int align_idx;
                int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                dpu_result_coord_t coord_var = {
                    .seq_nr = seq_number,
                    .seed_nr = sequence_idx,
                };

                if (seed_code < 0) {
                    continue;
                }

                seed = index_seed[seed_code];
                while (seed != NULL) {
                    if (seed_counter[seed_code].nb_seed < (int)seed->nb_nbr + total_nb_neighbour)
                        break;
                    total_nb_neighbour += seed->nb_nbr;
                    seed = seed->next;
                }
                align_idx = seed->offset + seed_counter[seed_code].nb_seed - total_nb_neighbour;

                code_neighbour(&ref_genome->data[sequence_start_idx + sequence_idx + SIZE_SEED], buf_code_neighbour);
                if (sequence_idx % 1000 == 0)
                    printf("\r\t%lli/%lli %i/%i         ", (unsigned long long)sequence_idx,
                        (unsigned long long)ref_genome->len_seq[seq_number] - SIZE_NEIGHBOUR_IN_BYTES - SIZE_SEED + 1, seq_number,
                        ref_genome->nb_seq);
                write_vmi(seed->num_dpu, align_idx, buf_code_neighbour, coord_var);

                seed_counter[seed_code].nb_seed++;
            }
            printf("\n");
        }
        free_vmis(nb_dpu);
        free(seed_counter);
        printf("\t\ttime: %lf s\n", my_clock() - write_in_memories_time);
    }

    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void index_free()
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
