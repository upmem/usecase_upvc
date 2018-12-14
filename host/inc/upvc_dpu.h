/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_DPU_H__
#define __UPVC_DPU_H__

#include <stdio.h>
#include <stdint.h>
#include "upvc.h"

#include "common.h"

/**
 * @brief DPU memory layout
 *
 * @var neighbour_idx        datas : table of the neighbours of the seed from the reference genome in the DPU
 * @var coordinate           datas : [sequence id, offset in sequence] of each seed in the DPU
 * @var neighbour_idx_len    Size of neighbour_idx, in bytes
 * @var coordinate_len       Size of coordinate, in bytes
 * @var neighbour_read       datas : table of neighbours
 * @var offset               datas : address of the first neighbour
 * @var count                datas : number of neighbour for each read
 * @var num                  datas : reads id
 */
typedef struct {
        int8_t *neighbour_idx;
        long *coordinate;
        unsigned int neighbour_idx_len;
        unsigned int coordinate_len;

        int8_t *neighbour_read;
        int *offset;
        int *count;
        int *num;
} mem_dpu_t;

mem_dpu_t *get_mem_dpu(unsigned int dpu_number);
dpu_result_out_t *get_mem_dpu_res(unsigned int dpu_number);

/**
 * @brief Allocate structure to store information of DPU memory.
 *
 * @param reads_info  Information on the size of the seed and the neighbour.
 * @param nb_dpu      Number of dpu used to compute.
 */
void malloc_dpu(reads_info_t *reads_info, int nb_dpu);
/**
 * @brief Allocate sub-Structure to store information of DPU memory.
 *
 * @param num_dpu     Number of the dpu for which to allocate the sub-structure.
 * @param nb_index    Number of index for this DPU.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
void malloc_neighbour_idx (int num_dpu, int nb_index, reads_info_t *reads_info);
/**
 * @brief Allocate structure to store results information of DPUs.
 *
 * @param nb_dpu  Number of dpu used to compute.
 */
void malloc_dpu_res(int nb_dpu);

/**
 * @brief Free all the structure allocated to store information of the DPU memory.
 *
 * @param nb_dpu  Number of DPUs.
 */
void free_dpu(int nb_dpu);
/**
 * @brief Free the structure used to store results information of DPUs.
 */
void free_dpu_res();

/**
 * @brief Write information of the DPU memory.
 */
void write_neighbours_and_coordinates(int d, int nb_nbrs, long *nbrs, reads_info_t *reads_info);
void write_neighbour_idx  (int num_dpu, int index_idx, int8_t *values, reads_info_t *reads_info);
void write_coordinate     (int num_dpu, int read_idx, long value);
void write_neighbour_read (int num_dpu, int read_idx, int8_t *values, reads_info_t *reads_info);
void write_count          (int num_dpu, int read_idx, int value);
void write_offset         (int num_dpu, int read_idx, int value);
void write_num            (int num_dpu, int read_idx, int value);

/**
 * @brief Read information of the DPU memory.
 */
int read_out_num          (int num_dpu, int align_idx);
int read_out_score        (int num_dpu, int align_idx);
dpu_result_coord_t read_out_coord       (int num_dpu, int align_idx);

/**
 * @brief Print information of the DPU memory.
 */
void print_neighbour_idx(int d, int offs, int nb_nbr, FILE *out, reads_info_t *reads_info);
void print_coordinates(int d, int offs, int l, FILE *out);

#endif /* __UPVC_DPU_H__ */
