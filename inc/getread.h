#ifndef __GETREAD_H__
#define __GETREAD_H__

#include <stdint.h>
#include "upvc.h"

/*
 * Fill "reads_buffer" with pairs of read.
 */
int  get_reads(FILE *fpe1, FILE *fpe2, int8_t *reads_buffer, times_ctx_t *times_ctx, reads_info_t *reads_info);
/*
 * Get the size of one read.
 */
int  get_read_size(char *filename_prefix);

#endif /* __GETREAD_H__ */
