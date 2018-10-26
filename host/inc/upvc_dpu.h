/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_DPU_H__
#define __UPVC_DPU_H__

#include <stdint.h>
#include "upvc.h"
#include "vmi.h"

/**
 * @brief Arguments needed for the thread that will be offload onto a DPU.
 */
typedef struct align_on_dpu_arg {
        reads_info_t reads_info;
        int numdpu;
} align_on_dpu_arg_t;

/**
 * @DPUs main function.
 */
void *align_on_dpu(void *arg);

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
 * @brief Free all the structure allocated to store information of the DPU memory.
 *
 * @param nb_dpu  Number of DPUs.
 */
void free_dpu(int nb_dpu);

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
long read_out_coord       (int num_dpu, int align_idx);

/**
 * @brief Print information of the DPU memory.
 */
void print_neighbour_idx(int d, int offs, int nb_nbr, FILE *out, reads_info_t *reads_info);
void print_coordinates(int d, int offs, int l, FILE *out);

/**
 * @brief Save the DPUs memory into files.
 */
void dump_mdpu_images_into_mram_files(vmi_t *vmis, unsigned int *nbr_counts, unsigned int nb_dpu, reads_info_t *reads_info);

/**
 * @brief Print a result for a DPU at a index.
 */
void write_result(int d, int k, unsigned int num, long coords, unsigned int score);

#endif /* __UPVC_DPU_H__ */
