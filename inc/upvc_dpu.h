#ifndef __UPVC_DPU_H__
#define __UPVC_DPU_H__

#include <stdint.h>
#include "upvc.h"

/*
 * Arguments needed for the thread that will be offload onto a DPU
 */
typedef struct align_on_dpu_arg {
        reads_info_t reads_info;
        int numdpu;
} align_on_dpu_arg_t;

/*
 * DPUs main function
 */
void *align_on_dpu(void *arg);

/*
 * Allocate DPU memory (TODO: to be erase once offloading on DPU is implemented)
 */
void malloc_dpu(reads_info_t *reads_info);
void malloc_neighbour_idx (int num_dpu, int nb_index, reads_info_t *reads_info);
/*
 * Free DPU memory
 */
void free_dpu();

/* Memory transfer with DPUs memories */
void write_neighbour_idx  (int num_dpu, int index_idx, int8_t *values, reads_info_t *reads_info);
void write_coordinate     (int num_dpu, int read_idx, long value);
void write_neighbour_read (int num_dpu, int read_idx, int8_t *values, reads_info_t *reads_info);
void write_count          (int num_dpu, int read_idx, int value);
void write_offset         (int num_dpu, int read_idx, int value);
void write_num            (int num_dpu, int read_idx, int value);
int read_out_num          (int num_dpu, int align_idx);
int read_out_score        (int num_dpu, int align_idx);
long read_out_coord       (int num_dpu, int align_idx);

#endif /* __UPVC_DPU_H__ */
