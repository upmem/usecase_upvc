/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __SIMU_BACKEND_H__
#define __SIMU_BACKEND_H__

#include <semaphore.h>

void run_dpu_simulation(unsigned int dpu_offset, unsigned int pass_id, sem_t *dispatch_free_sem, sem_t *acc_wait_sem,
    sem_t *exec_to_acc_sem, sem_t *dispatch_to_exec_sem);

void init_backend_simulation(unsigned int *nb_dpus_per_run);

void free_backend_simulation();

void load_mram_simulation(unsigned int dpu_offset, int delta_neighbour);

void wait_dpu_simulation();

#endif /* __SIMU_BACKEND_H__ */
