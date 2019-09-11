/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __DISPATCH_H__
#define __DISPATCH_H__

#include <stdint.h>

/**
 * @brief List of reads dispatched to a DPU.
 *
 * @var nb_reads  The number of requests.
 * @var reads     A table of nb_reads requests. Since the read size is not fixed, the table is a raw byte stream.
 */
typedef struct {
    uint32_t nb_reads;
    int8_t *reads_area;
} dispatch_request_t;

#include "dpus_mgmt.h"
#include "index.h"
#include "upvc.h"
#include "vmi.h"

#include "backends_functions.h"

/**
 * @brief Dispatch reads amongst the list of DPUs, to create requests.
 *
 * @param index_seed       List of the seed of the reference genome.
 * @param read_buffer      Buffer containning the reads to dispatch.
 * @param nb_read          Number of reads to dispatch.
 * @param dispatch_requests  Table of dispatch_request_t to store the dispatching between the DPUs.
 * @param times_ctx        Times information for the whole application.
 * @param reads_info       Information on the size of the seed and the neighbour.
 * @param backends_functions  Functions to be use by the dispatcher depending on the backend.
 *
 * @return The dispatcher result.
 */
void dispatch_read(index_seed_t **index_seed, int8_t *read_buffer, int nb_read, dispatch_request_t *dispatch_requests,
    times_ctx_t *times_ctx, reads_info_t *reads_info, backends_functions_t *backends_functions);

/**
 * @brief Create a table of requests to be filled with dispatch_read
 *
 * @param nb_dpu      Number of DPUs on which to dispatch the reads.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
dispatch_request_t *dispatch_create(unsigned int nb_dpu, reads_info_t *reads_info);

/**
 * @brief Frees the requests produced by dispatch_read.
 *
 * @param dispatch  The free structure.
 * @param nb_dpu    The number of DPUs to use to compute.
 */
void dispatch_free(dispatch_request_t *dispatch, unsigned int nb_dpu);

#endif /* __DISPATCH_H__ */
