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
    void (*run_dpu)(unsigned int, unsigned int, sem_t *, sem_t *, sem_t *, sem_t *);
    void (*init_backend)(unsigned int *);
    void (*free_backend)(void);
    void (*load_mram)(unsigned int, int);
    void (*wait_dpu)(void);
} backends_functions_t;

#endif /* __BACKENDS_FUNCTIONS_H__ */
