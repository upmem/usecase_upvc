#ifndef __DISPATCH_H__
#define __DISPATCH_H__

#include <stdint.h>
#include "index.h"
#include "upvc.h"

/*
 * Dispatch the seed from the reference genome onto the DPUs
 */
void dispatch_read(index_seed_t **index_seed, int8_t *read_buffer, int nb_read, times_ctx_t *times_ctx, reads_info_t *reads_info);

#endif /* __DISPATCH_H__ */
