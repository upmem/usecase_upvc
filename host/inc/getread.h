/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __GETREAD_H__
#define __GETREAD_H__

#include "upvc.h"
#include <stdint.h>
#include <stdio.h>

/**
 * @brief Fill "reads_buffer" with pairs of read.
 *
 * @param fpe1          First file of Pair-end read.
 * @param fpe2          Second file of Pair-end read.
 * @param reads_buffer  Output buffer to save pairs of read.
 * @param times_ctx     Times information for the while application.
 * @param reads_info    Information on the size of the seed and the neighbour.
 *
 * @return The number of read written in "reads_buffer".
 */
int get_reads(FILE *fpe1, FILE *fpe2, int8_t *reads_buffer, times_ctx_t *times_ctx, reads_info_t *reads_info);

/**
 * @brief Get the size of one read.
 *
 * @param input_pe1_file  File to parse to find the size of one read.
 *
 * @return The size of one read.
 */
int get_read_size(char *input_pe1_file);

#endif /* __GETREAD_H__ */
