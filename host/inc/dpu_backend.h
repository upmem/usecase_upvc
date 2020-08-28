/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __DPU_BACKEND_H__
#define __DPU_BACKEND_H__

#include <semaphore.h>

void get_dpu_id(unsigned int numdpu, unsigned int *, unsigned int *, unsigned int *);

void run_on_dpu(
    unsigned int dpu_offset, unsigned int rank_id, unsigned int pass_id, sem_t *dispatch_free_sem, sem_t *acc_wait_sem);

void init_backend_dpu(unsigned int *nb_dpus_per_run, unsigned int *nb_ranks_per_run);

void free_backend_dpu();

void load_mram_dpu(unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour);

int get_rank_numa_node(int rank_id);

#endif /* __DPU_BACKEND_H__ */
