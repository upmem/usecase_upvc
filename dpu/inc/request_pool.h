/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __REQUEST_POOL_H__
#define __REQUEST_POOL_H__

#include <mram.h>

#include "common.h"
#include "dout.h"

__mram_ptr void *get_request_addr();
size_t get_request_size();

/**
 * @brief Gets the next read from the request pool, if any.
 *
 * @param request_buffer  Output Request filled with the next request and corresponding neighbour.
 * @param stats           To update statistical reports.
 *
 * @return True if a new request was fetched, false if the FIFO is empty.
 */
bool request_pool_next(dpu_request_t *request_buffer, dpu_tasklet_stats_t *stats);

/**
 * @brief Initializes the request pool.
 */
void request_pool_init();

#endif /* __REQUEST_POOL_H__ */
