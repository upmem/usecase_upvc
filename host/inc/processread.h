/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __PROCESSREAD_H__
#define __PROCESSREAD_H__

#include <stdio.h>

typedef struct {
    int type;
    int ix;
    int jx;
} backtrack_t;


#define CODE_MATCH 0
#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13
#define CODE_ERR 14

void process_read(FILE *fpe1, FILE *fpe2, int round, unsigned int pass_id);

void process_read_init();

void process_read_free();

#endif /* __PROCESSREAD_H__ */
