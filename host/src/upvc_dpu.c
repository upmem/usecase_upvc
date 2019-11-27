/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upvc.h"
#include "upvc_dpu.h"

#include "common.h"

static mem_dpu_t *MDPU;
static dpu_result_out_t **MDPU_res;

mem_dpu_t *get_mem_dpu(unsigned int dpu_number) { return &MDPU[dpu_number]; }

dpu_result_out_t *get_mem_dpu_res(unsigned int dpu_number) { return MDPU_res[dpu_number]; }

void malloc_dpu_res(int nb_dpu)
{
    MDPU_res = (dpu_result_out_t **)malloc(sizeof(dpu_result_out_t *) * nb_dpu);
    for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
        MDPU_res[num_dpu] = (dpu_result_out_t *)malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
    }
}

void malloc_dpu(int nb_dpu)
{
    malloc_dpu_res(nb_dpu);
    MDPU = (mem_dpu_t *)calloc(nb_dpu, sizeof(mem_dpu_t));
    assert(MDPU != NULL);

    for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
        MDPU[num_dpu].neighbour_read = (int8_t *)malloc(sizeof(int8_t) * SIZE_NEIGHBOUR_IN_BYTES * MAX_DPU_REQUEST);
        MDPU[num_dpu].count = (int *)malloc(sizeof(int) * MAX_DPU_REQUEST);
        MDPU[num_dpu].offset = (int *)malloc(sizeof(int) * MAX_DPU_REQUEST);
        MDPU[num_dpu].num = (int *)malloc(sizeof(int) * MAX_DPU_REQUEST);
        MDPU[num_dpu].neighbour_idx = (int8_t *)malloc(MRAM_SIZE);
        assert(MDPU[num_dpu].neighbour_idx != NULL);
        assert(MDPU[num_dpu].neighbour_read != NULL);
        assert(MDPU[num_dpu].count != NULL);
        assert(MDPU[num_dpu].offset != NULL);
        assert(MDPU[num_dpu].num != NULL);
    }
}

void free_dpu_res(int nb_dpu)
{
    for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
        free(MDPU_res[num_dpu]);
    }
    free(MDPU_res);
}

void free_dpu(int nb_dpu)
{
    for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
        free(MDPU[num_dpu].neighbour_idx);
        free(MDPU[num_dpu].neighbour_read);
        free(MDPU[num_dpu].count);
        free(MDPU[num_dpu].offset);
        free(MDPU[num_dpu].num);
    }
    free_dpu_res(nb_dpu);
    free(MDPU);
}

void write_neighbours_and_coordinates(int numdpu, int index_idx, int8_t *nbrs, dpu_result_coord_t coord)
{
    mem_dpu_t *mdpu = &(MDPU[numdpu]);

    unsigned int size_neighbour_in_bytes = SIZE_NEIGHBOUR_IN_BYTES;
    unsigned int aligned_nbr_len = ALIGN_DPU(size_neighbour_in_bytes);
    unsigned int aligned_coord_nbr_len = sizeof(dpu_result_coord_t) + aligned_nbr_len;

    *(dpu_result_coord_t *)(mdpu->neighbour_idx + (index_idx * aligned_coord_nbr_len)) = coord;
    (void)memcpy(
        mdpu->neighbour_idx + (index_idx * aligned_coord_nbr_len) + sizeof(dpu_result_coord_t), nbrs, (size_t)aligned_nbr_len);
}

void write_neighbour_read(int num_dpu, int read_idx, int8_t *values)
{
    int size_neighbour = SIZE_NEIGHBOUR_IN_BYTES;
    for (int i = 0; i < size_neighbour; i++) {
        MDPU[num_dpu].neighbour_read[read_idx * size_neighbour + i] = values[i];
    }
}

void write_count(int num_dpu, int read_idx, int value) { MDPU[num_dpu].count[read_idx] = value; }
void write_offset(int num_dpu, int read_idx, int value) { MDPU[num_dpu].offset[read_idx] = value; }
void write_num(int num_dpu, int read_idx, int value) { MDPU[num_dpu].num[read_idx] = value; }

int read_out_num(int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].num; }
int read_out_score(int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].score; }
dpu_result_coord_t read_out_coord(int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].coord; }

static void print_byte_to_sym(uint8_t byte, unsigned int offset, FILE *out)
{
    uint8_t as_byte = (uint8_t)((byte >> offset) & 0x3);
    char as_char;
    switch (as_byte) {
    case CODE_A:
        as_char = 'A';
        break;
    case CODE_C:
        as_char = 'C';
        break;
    case CODE_G:
        as_char = 'G';
        break;
    default:
        as_char = 'T';
    }
    fprintf(out, "%c", as_char);
}

void print_neighbour_idx(int d, int offs, int nb_nbr, FILE *out)
{
    unsigned int size_neighbour_in_bytes = SIZE_NEIGHBOUR_IN_BYTES;
    unsigned int aligned_nbr_len = ALIGN_DPU(size_neighbour_in_bytes);
    unsigned int aligned_coord_nbr_len = sizeof(dpu_result_coord_t) + aligned_nbr_len;
    for (int each_nbr = 0; each_nbr < nb_nbr; each_nbr++) {
        fprintf(out, "\t");
        unsigned int i;
        for (i = 0; i < size_neighbour_in_bytes; i++) {
            uint8_t this_byte
                = (uint8_t)MDPU[d].neighbour_idx[(offs + each_nbr) * aligned_coord_nbr_len + sizeof(dpu_result_coord_t) + i];
            print_byte_to_sym(this_byte, 0, out);
            print_byte_to_sym(this_byte, 2, out);
            print_byte_to_sym(this_byte, 4, out);
            print_byte_to_sym(this_byte, 6, out);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

void print_coordinates(int d, int offs, int l, FILE *out)
{
    unsigned int size_neighbour_in_bytes = SIZE_NEIGHBOUR_IN_BYTES;
    unsigned int aligned_nbr_len = ALIGN_DPU(size_neighbour_in_bytes);
    unsigned int aligned_coord_nbr_len = sizeof(dpu_result_coord_t) + aligned_nbr_len;
    for (int i = 0; i < l; i++) {
        fprintf(out, " %llu",
            (unsigned long long)((dpu_result_coord_t *)(&MDPU[d].neighbour_idx[(offs + l) * aligned_coord_nbr_len]))->coord);
    }
    fprintf(out, "\n");
}
