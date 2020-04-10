/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __SIMU_BACKEND_H__
#define __SIMU_BACKEND_H__

#include <semaphore.h>

#include "dispatch.h"
#include "index.h"

/**
 * @brief Write seed information in the structure representing the DPU list of requests.
 */
void add_seed_to_simulation_requests(
    dispatch_request_t *requests, int num_read, int nb_read_written, index_seed_t *seed, int8_t *nbr);

/**
 * @brief Compute one pass in simulation mode.
 */
void run_dpu_simulation(unsigned int dpu_offset, unsigned int rank_id, unsigned int pass_id, int delta_neighbour,
    sem_t *dispatch_free_sem, sem_t *acc_wait_sem);

/**
 * @brief Index the reference genome.
 */
void init_backend_simulation();

/**
 * @brief Free structure used by simulation.
 */
void free_backend_simulation(unsigned int nb_dpu);

/**
 * @brief Load mram in structure that will represente the DPUs.
 */
void load_mram_simulation(unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour);

#endif /* __SIMU_BACKEND_H__ */
