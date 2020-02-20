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
 * @var nb_reads      The number of reads in the request pool.
 * @var rdidx         Index of the first unread read in the request pool.
 * @var cur_read      Address of the first read to be processed in MRAM.
 */
typedef struct {
    unsigned int nb_reads;
    unsigned int rdidx;
    uintptr_t cur_read;
} request_pool_t;

__mram_noinit nb_request_t DPU_NB_REQUEST_VAR;

__mram_noinit dpu_request_t DPU_REQUEST_VAR[MAX_DPU_REQUEST];

/**
 * @brief Common request pool, shared by every tasklet.
 */
static request_pool_t request_pool;
MUTEX_INIT(request_pool_mutex);

void request_pool_init()
{
    request_pool.nb_reads = (unsigned int)DPU_NB_REQUEST_VAR;
    request_pool.rdidx = 0;
    request_pool.cur_read = (uintptr_t)DPU_REQUEST_VAR;
}

bool request_pool_next(dpu_request_t *request, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    mutex_lock(request_pool_mutex);
    if (request_pool.rdidx == request_pool.nb_reads) {
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
