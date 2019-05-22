/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __DPU_BACKEND_H__
#define __DPU_BACKEND_H__

#include <semaphore.h>

#include "genome.h"
#include "vmi.h"
#include "upvc.h"
#include "index.h"
#include "dispatch.h"
#include "backends_functions.h"

/**
 * @brief Add a seed to a request.
 */
void add_seed_to_dpu_requests(dispatch_request_t *requests,
                              int num_read,
                              int nb_read_written,
                              index_seed_t *seed,
                              int8_t *nbr,
                              reads_info_t *reads_info);

/**
 * @brief Compute one pass on DPUs.
 */
void run_on_dpu(dispatch_request_t *dispatch,
                devices_t *devices,
                unsigned int dpu_offset,
                unsigned int rank_id,
                unsigned int round,
                unsigned int nb_pass,
                sem_t *dispatch_free_sem,
                sem_t *acc_wait_sem,
                times_ctx_t *times_ctx,
                reads_info_t *reads_info);

/**
 * @brief Index the reference genome and alloc physical DPUs
 */
void init_backend_dpu(unsigned int *nb_rank,
                      devices_t **devices,
                      unsigned int nb_dpu_per_run,
                      const char *dpu_binary,
                      index_seed_t ***index_seed,
                      reads_info_t *reads_info);

/**
 * @brief Free DPUs.
 */
void free_backend_dpu(devices_t *devices, unsigned int nb_dpu);

/**
 * @brief load mram into DPUs for one run
 */
void load_mram_dpu(unsigned int dpu_offset,
                   unsigned int rank_id,
                   devices_t *devices,
                   reads_info_t *reads_info,
                   times_ctx_t *times_ctx);

#endif /* __DPU_BACKEND_H__ */
