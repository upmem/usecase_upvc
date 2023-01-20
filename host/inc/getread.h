/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __GETREAD_H__
#define __GETREAD_H__

#include <stdint.h>
#include <stdio.h>

int get_reads_in_buffer(unsigned int pass_id);

int8_t *get_reads_buffer(unsigned int pass_id);
float *get_reads_quality_buffer(unsigned int pass_id);
char *get_reads_name_buffer(unsigned int pass_id);

void get_reads(FILE *fpe1, FILE *fpe2, unsigned int pass_id);

int get_input_info(FILE *f, size_t *read_size, size_t *nb_read);

#endif /* __GETREAD_H__ */
