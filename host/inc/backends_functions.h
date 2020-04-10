/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __BACKENDS_FUNCTIONS_H__
#define __BACKENDS_FUNCTIONS_H__

#include <semaphore.h>

#include "dispatch.h"
#include "index.h"

/**
 * @brief Specific function whether the application run on simulation or on DPU (fsim of fpga).
 */
typedef struct backends_functions_struct {
    void (*run_dpu)(unsigned int, unsigned int, unsigned int, int, sem_t *, sem_t *);
    void (*add_seed_to_requests)(dispatch_request_t *, int, int, index_seed_t *, int8_t *);

    void (*init_backend)(void);
    void (*free_backend)(unsigned int);
    void (*load_mram)(unsigned int, unsigned int, int);
} backends_functions_t;

extern backends_functions_t backends_functions;

#endif /* __BACKENDS_FUNCTIONS_H__ */
