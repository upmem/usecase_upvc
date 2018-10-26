/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __DISPATCH_H__
#define __DISPATCH_H__

#include <stdbool.h>
#include <stdint.h>
#include "index.h"
#include "upvc.h"

/**
 * @brief Structure representing one request to a DPU, resulting from the dispatching of reads.
 *
 * Such a request basically contains all the information for one read.
 *
 * @var offset  The 1st neighbour address.
 * @var count   The number of neighbours.
 * @var num     A reference number to the original request.
 * @var nbr     The neighbours part of this read, actual size of this field is variable.
 */
typedef struct {
        uint32_t offset;
        uint32_t count;
        uint32_t num;
        uint32_t nbr;
} dispatch_read_t;

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

/**
 * @brief A list of requests distributed amongst DPUs.
 *
 * A table of requests, indexed by target DPU.
 */
typedef dispatch_request_t *dispatch_t;

/**
 * @brief Dispatch reads amongst the list of DPUs, to create requests.
 *
 * @param index_seed       List of the seed of the reference genome.
 * @param read_buffer      Buffer containning the reads to dispatch.
 * @param nb_read          Number of reads to dispatch.
 * @param nb_dpu           Number of DPUs on which to dispatch the reads.
 * @param times_ctx        Times information for the whole application.
 * @param reads_info       Information on the size of the seed and the neighbour.
 * @param simulation_mode  Indicate if we have DPUs (fpga of hsim) to compute or only the host.
 *
 * @return The dispatcher result.
 */
dispatch_t dispatch_read(index_seed_t **index_seed,
                         int8_t *read_buffer,
                         int nb_read,
                         int nb_dpu,
                         times_ctx_t *times_ctx,
                         reads_info_t *reads_info,
                         bool simulation_mode);

/**
 * @brief Frees the requests produced by dispatch_read.
 *
 * @param dispatch  The free structure.
 * @param nb_dpu    The number of DPUs to use to compute.
 */
void dispatch_free(dispatch_t dispatch, unsigned int nb_dpu);

/**
 * @brief Get the exact size, in bytes, of a read with this execution context.
 *
 * @param reads_info  Information on the size of the seed and the neighbour.
 *
 * @return The size of a read.
 */
unsigned int dispatch_read_len(reads_info_t *reads_info);

#endif /* __DISPATCH_H__ */
