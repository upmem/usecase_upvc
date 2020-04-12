/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __DPU_BACKEND_H__
#define __DPU_BACKEND_H__

#include <semaphore.h>

#include "dispatch.h"
#include "index.h"

void run_on_dpu(unsigned int dpu_offset, unsigned int rank_id, unsigned int pass_id, int delta_neighbour,
    sem_t *dispatch_free_sem, sem_t *acc_wait_sem);

void init_backend_dpu();

void free_backend_dpu();

void load_mram_dpu(unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour);

unsigned int get_nb_ranks_per_run_dpu();

#endif /* __DPU_BACKEND_H__ */
