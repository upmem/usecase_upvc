/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __SIMU_BACKEND_H__
#define __SIMU_BACKEND_H__

#include <semaphore.h>

#include "backends_functions.h"
#include "dispatch.h"
#include "genome.h"
#include "index.h"
#include "upvc.h"
#include "vmi.h"

/**
 * @brief Write seed information in the structure representing the DPU list of requests.
 */
void add_seed_to_simulation_requests(
    dispatch_request_t *requests, int num_read, int nb_read_written, index_seed_t *seed, int8_t *nbr);

/**
 * @brief Compute one pass in simulation mode.
 */
void run_dpu_simulation(dispatch_request_t *dispatch, devices_t *devices, unsigned int dpu_offset, unsigned int rank_id,
    int delta_neighbour, sem_t *dispatch_free_sem, sem_t *acc_wait_sem, times_ctx_t *times_ctx);

/**
 * @brief Index the reference genome.
 */
void init_backend_simulation(
    unsigned int *nb_rank, devices_t **devices, unsigned int nb_dpu_per_run, index_seed_t ***index_seed);

/**
 * @brief Free structure used by simulation.
 */
void free_backend_simulation(devices_t *devices, unsigned int nb_dpu);

/**
 * @brief Load mram in structure that will represente the DPUs.
 */
void load_mram_simulation(
    unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour, devices_t *devices, times_ctx_t *times_ctx);

#endif /* __SIMU_BACKEND_H__ */
