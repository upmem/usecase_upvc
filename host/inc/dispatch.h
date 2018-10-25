#ifndef __DISPATCH_H__
#define __DISPATCH_H__

#include <stdbool.h>
#include <stdint.h>
#include "index.h"
#include "upvc.h"

/**
 * @brief Structure representing one request to a DPU, resulting from the dispatching of reads
 *
 * Such a request basically contains all the information for one read.
 *
 * @var offset  the 1st neighbor address
 * @var count   the number of neighbors
 * @var num     a reference number to the original request
 * @var nbr     the neighbors part of this read, actual size of this field is variable
 */
typedef struct {
    uint32_t offset;
    uint32_t count;
    uint32_t num;
    uint32_t nbr;
} dispatch_read_t;

/**
 * @brief List of reads dispatched to a DPU
 * @var nr_reads  the number of requests
 * @var reads     a table of nr_reads requests. Since the read size is not fixed, the table is a raw byte stream.
 */
typedef struct {
    uint32_t nr_reads;
    int8_t *reads_area;
} dispatch_request_t;

/**
 * @brief A list of requests distributed amongst DPUs.
 *
 * A table of requests, indexed by target DPU.
 */
typedef dispatch_request_t *dispatch_t;

/*
 * Dispatch the seed from the reference genome onto the DPUs
 */
dispatch_t dispatch_read(index_seed_t **index_seed,
                         int8_t *read_buffer,
                         int nb_read,
                         int nb_dpu,
                         times_ctx_t *times_ctx,
                         reads_info_t *reads_info,
                         bool simulation_mode);

/**
 * @brief Frees the requests produced by dispatch_read
 * @param dispatch  the free structure
 * @param nb_dpu is the number of dpu to use to compute
 */
void dispatch_free(dispatch_t dispatch, unsigned int nb_dpu);

/**
 * @param reads_info are informations about seed and neighbour size.
 * @return The exact size, in bytes, of a read within this execution context.
 */
unsigned int dispatch_read_len(reads_info_t *reads_info);

#endif /* __DISPATCH_H__ */
