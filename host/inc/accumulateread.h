/**
 * Copyright 2020 - Dominique Lavenier & UPMEM
 */

#ifndef __ACCUMULATEREAD_H__
#define __ACCUMULATEREAD_H__

#include "common.h"
#include "stdbool.h"

typedef struct {
    nb_result_t nb_res;
    dpu_result_out_t *results;
} acc_results_t;

acc_results_t *accumulate_get_buffer(unsigned int dpu_id, unsigned int pass_id);
acc_results_t accumulate_get_result(unsigned int pass_id);


void accumulate_disable_dpus(void);
void accumulate_read(unsigned int pass_id, unsigned int max_nb_pass);
uint32_t accumulate_valid_run();

void accumulate_init();
void accumulate_free();

#endif /* __ACCUMULATEREAD_H__ */
