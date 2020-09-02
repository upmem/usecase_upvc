/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <string.h>

#include <alloc.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>

#include "debug.h"
#include "dout.h"
#include "macro_utils.h"
#include "request_pool.h"
#include "stats.h"

#include "common.h"

/**
 * @brief Common structure to consume requests.
 *
 * Requests belong to a FIFO, from which each tasklet picks the reads. The FIFO is protected by a critical
 * section.
 *
 * @var rdidx         Index of the first unread read in the request pool.
 * @var cur_read      Address of the first read to be processed in MRAM.
 */
typedef struct {
    uint32_t rdidx;
    uintptr_t cur_read;
} request_pool_t;

__host nb_request_t DPU_NB_REQUEST_VAR0;
__host nb_request_t DPU_NB_REQUEST_VAR1;
__host nb_request_t DPU_NB_REQUEST_VAR2;
__host nb_request_t DPU_NB_REQUEST_VAR3;
nb_request_t *DPU_NB_REQUEST_VAR;

__mram_noinit dpu_request_t DPU_REQUEST_VAR0[MAX_DPU_REQUEST];
__mram_noinit dpu_request_t DPU_REQUEST_VAR1[MAX_DPU_REQUEST];
__mram_noinit dpu_request_t DPU_REQUEST_VAR2[MAX_DPU_REQUEST];
__mram_noinit dpu_request_t DPU_REQUEST_VAR3[MAX_DPU_REQUEST];

__mram_ptr dpu_request_t *DPU_REQUEST_VAR;
extern uint32_t dpu_result_var_idx;

/**
 * @brief Common request pool, shared by every tasklet.
 */
static request_pool_t request_pool;
MUTEX_INIT(request_pool_mutex);

void request_pool_init()
{
    if (dpu_result_var_idx % 4 == 0) {
        DPU_NB_REQUEST_VAR = &DPU_NB_REQUEST_VAR0;
        DPU_REQUEST_VAR = DPU_REQUEST_VAR0;
    } else if (dpu_result_var_idx % 4 == 1) {
        DPU_NB_REQUEST_VAR = &DPU_NB_REQUEST_VAR1;
        DPU_REQUEST_VAR = DPU_REQUEST_VAR1;
    } else if (dpu_result_var_idx % 4 == 2) {
        DPU_NB_REQUEST_VAR = &DPU_NB_REQUEST_VAR2;
        DPU_REQUEST_VAR = DPU_REQUEST_VAR2;
    } else if (dpu_result_var_idx % 4 == 3) {
        DPU_NB_REQUEST_VAR = &DPU_NB_REQUEST_VAR3;
        DPU_REQUEST_VAR = DPU_REQUEST_VAR3;
    }

    request_pool.rdidx = 0;
    request_pool.cur_read = (uintptr_t)DPU_REQUEST_VAR;
}

bool request_pool_next(dpu_request_t *request, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    mutex_lock(request_pool_mutex);
    if (request_pool.rdidx == *DPU_NB_REQUEST_VAR) {
        mutex_unlock(request_pool_mutex);
        return false;
    }

    /* Fetch next request into cache */
    STATS_INCR_LOAD(stats, sizeof(dpu_request_t));
    STATS_INCR_LOAD_DATA(stats, sizeof(dpu_request_t));
    mram_read((__mram_ptr void *)request_pool.cur_read, (void *)request, sizeof(*request));

    /* Point to next request */
    request_pool.rdidx++;
    request_pool.cur_read += sizeof(dpu_request_t);
    mutex_unlock(request_pool_mutex);

    return true;
}
