/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __BACKENDS_FUNCTIONS_H__
#define __BACKENDS_FUNCTIONS_H__

#include <semaphore.h>

/**
 * @brief Specific function whether the application run on simulation or on DPU (fsim of fpga).
 */
typedef struct backends_functions_struct {
    void (*run_dpu)(dispatch_request_t *, devices_t *, unsigned int, unsigned int, unsigned int, sem_t *, sem_t *, times_ctx_t *,
        reads_info_t *);
    void (*add_seed_to_requests)(dispatch_request_t *, int, int, index_seed_t *, int8_t *, reads_info_t *);

    void (*init_backend)(unsigned int *, devices_t **, unsigned int, const char *, index_seed_t ***, reads_info_t *);
    void (*free_backend)(devices_t *, unsigned int);
    void (*load_mram)(unsigned int, unsigned int, devices_t *, reads_info_t *, times_ctx_t *);
} backends_functions_t;

#endif /* __BACKENDS_FUNCTIONS_H__ */
