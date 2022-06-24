/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __PROCESSREAD_H__
#define __PROCESSREAD_H__

#include <stdio.h>

void process_read(FILE *fpe1, FILE *fpe2, int round, unsigned int pass_id);

void process_read_init();

void process_read_free();

#endif /* __PROCESSREAD_H__ */
