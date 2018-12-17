/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "upvc_dpu.h"
#include "upvc.h"

#include "common.h"

static mem_dpu_t *MDPU;
static dpu_result_out_t **MDPU_res;

mem_dpu_t *get_mem_dpu(unsigned int dpu_number)
{
        return &MDPU[dpu_number];
}

dpu_result_out_t *get_mem_dpu_res(unsigned int dpu_number)
{
        return MDPU_res[dpu_number];
}

void malloc_dpu_res(int nb_dpu) {
        MDPU_res = (dpu_result_out_t **) malloc(sizeof(dpu_result_out_t*) * nb_dpu);
        for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
                MDPU_res[num_dpu] = (dpu_result_out_t *) malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
        }
}

void malloc_dpu(reads_info_t *reads_info, int nb_dpu)
{
        malloc_dpu_res(nb_dpu);
        MDPU = (mem_dpu_t *) malloc(sizeof(mem_dpu_t) * nb_dpu);

        for (int num_dpu = 0; num_dpu < nb_dpu; num_dpu++) {
                MDPU[num_dpu].neighbour_idx = NULL;
                MDPU[num_dpu].neighbour_idx_len = 0;
                MDPU[num_dpu].coordinate = NULL;
                MDPU[num_dpu].coordinate_len = 0;
                MDPU[num_dpu].neighbour_read =
                        (int8_t *) malloc(sizeof(int8_t) * reads_info->size_neighbour_in_bytes * MAX_DPU_REQUEST);
                MDPU[num_dpu].count          = (int *)    malloc(sizeof(int) * MAX_DPU_REQUEST);
                MDPU[num_dpu].offset         = (int *)    malloc(sizeof(int) * MAX_DPU_REQUEST);
                MDPU[num_dpu].num            = (int *)    malloc(sizeof(int) * MAX_DPU_REQUEST);
        }
}

void malloc_neighbour_idx (int num_dpu, int nb_index, reads_info_t *reads_info)
{
        MDPU[num_dpu].neighbour_idx_len = sizeof(int8_t) * nb_index * reads_info->size_neighbour_in_bytes;
        MDPU[num_dpu].neighbour_idx = (int8_t *) malloc(MDPU[num_dpu].neighbour_idx_len);
        MDPU[num_dpu].coordinate_len = sizeof(long) * nb_index;
        MDPU[num_dpu].coordinate = (long *) malloc(MDPU[num_dpu].coordinate_len);
}

void free_dpu_res()
{
        free(MDPU_res);
}

void free_dpu(int nb_dpu)
{
        for (int num_dpu = 0; num_dpu<nb_dpu; num_dpu++) {
                free (MDPU[num_dpu].neighbour_idx);
                free (MDPU[num_dpu].neighbour_read);
                free (MDPU[num_dpu].coordinate);
                free (MDPU[num_dpu].count);
                free (MDPU[num_dpu].offset);
                free (MDPU[num_dpu].num);
                free (MDPU_res[num_dpu]);
        }
        free_dpu_res();
        free(MDPU);
}

void write_neighbours_and_coordinates(int d, int nb_nbrs, long *nbrs, reads_info_t *reads_info)
{
        mem_dpu_t *mdpu = &(MDPU[d]);

        unsigned int size_neighbour_in_bytes = reads_info->size_neighbour_in_bytes;
        size_t long_aligned_nbr_len = (size_t) ((size_neighbour_in_bytes + 7) & (~7));
        unsigned int input_pattern_len = sizeof(long) + long_aligned_nbr_len;
        uint8_t *input_pattern_ptr = (uint8_t *) nbrs;

        for (int each_nbr = 0; each_nbr < nb_nbrs; each_nbr++) {
                long this_coords = ((long *) input_pattern_ptr)[0];
                uint8_t *this_nbr = input_pattern_ptr + sizeof(long);
                (void) memcpy(mdpu->neighbour_idx + (each_nbr * size_neighbour_in_bytes),
                              this_nbr,
                              (size_t) size_neighbour_in_bytes);
                mdpu->coordinate[each_nbr] = this_coords;
                input_pattern_ptr += input_pattern_len;
        }
}

void write_neighbour_read (int num_dpu, int read_idx, int8_t *values, reads_info_t *reads_info)
{
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        for (int i = 0; i < size_neighbour; i++) {
                MDPU[num_dpu].neighbour_read[read_idx * size_neighbour + i] = values[i];
        }
}

void write_neighbour_idx (int num_dpu, int index_idx, int8_t *values, reads_info_t *reads_info)
{
        int size_neighbour = reads_info->size_neighbour_in_bytes;
        for (int i = 0; i < size_neighbour; i++) {
                MDPU[num_dpu].neighbour_idx[index_idx * size_neighbour + i] = values[i];
        }
}

void write_coordinate  (int num_dpu, int read_idx, long value) { MDPU[num_dpu].coordinate[read_idx] = value; }
void write_count       (int num_dpu, int read_idx, int value)  { MDPU[num_dpu].count[read_idx]      = value; }
void write_offset      (int num_dpu, int read_idx, int value)  { MDPU[num_dpu].offset[read_idx]     = value; }
void write_num         (int num_dpu, int read_idx, int value)  { MDPU[num_dpu].num[read_idx]        = value; }

int read_out_num       (int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].num;   }
int read_out_score     (int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].score; }
dpu_result_coord_t read_out_coord    (int num_dpu, int align_idx) { return MDPU_res[num_dpu][align_idx].coord; }

static void print_byte_to_sym(uint8_t byte, unsigned int offset, FILE *out)
{
        uint8_t as_byte = (uint8_t) ((byte >> offset) & 0x3);
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

void print_neighbour_idx(int d, int offs, int nb_nbr, FILE *out, reads_info_t *reads_info)
{
        unsigned int size_neighbour_in_bytes = reads_info->size_neighbour_in_bytes;
        for (int each_nbr = 0; each_nbr < nb_nbr; each_nbr++) {
                fprintf(out, "\t");
                unsigned int i;
                for (i = 0; i < size_neighbour_in_bytes; i++) {
                        uint8_t this_byte = (uint8_t) MDPU[d].neighbour_idx[(offs + each_nbr) * size_neighbour_in_bytes + i];
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
        for (int i = 0; i < l; i++) {
                fprintf(out, " %lu", MDPU[d].coordinate[offs + i]);
        }
        fprintf(out, "\n");
}
