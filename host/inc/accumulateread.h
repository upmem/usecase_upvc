/**
 * Copyright 2020 - Dominique Lavenier & UPMEM
 */

#ifndef __ACCUMULATEREAD_H__
#define __ACCUMULATEREAD_H__

#include "common.h"

void accumulate_read(unsigned int pass_id, unsigned int dpu_offset);
void accumulate_get_result(unsigned int pass_id, unsigned int *nb_res, dpu_result_out_t **results);

void accumulate_init();
void accumulate_free();

#endif /* __ACCUMULATEREAD_H__ */
