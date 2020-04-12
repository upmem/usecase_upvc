/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __DISPATCH_H__
#define __DISPATCH_H__

#include <stdint.h>

#include "common.h"

/**
 * @brief List of reads dispatched to a DPU.
 *
 * @var nb_reads  The number of requests.
 * @var reads     A table of nb_reads requests. Since the read size is not fixed, the table is a raw byte stream.
 */
typedef struct {
    nb_request_t nb_reads;
    dpu_request_t *dpu_requests;
} dispatch_request_t;

dispatch_request_t *dispatch_get(unsigned int dpu_id, unsigned int pass_id);

void dispatch_read(unsigned int pass_id);

void dispatch_init(unsigned int nb_dpu);
void dispatch_free(unsigned int nb_dpu);

#endif /* __DISPATCH_H__ */
