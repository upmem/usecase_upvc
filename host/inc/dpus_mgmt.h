/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DPUS_H__
#define __INTEGRATION_DPUS_H__

#include "common.h"

unsigned int dpu_get_nb_ranks_per_run();

unsigned int dpu_try_alloc_for();

void dpu_try_write_mram(unsigned int rank_id, unsigned int dpu_offset, delta_info_t delta_neighbour);

void dpu_try_free();

void dpu_try_run(unsigned int rank_id);

void dpu_try_write_dispatch_into_mram(unsigned int rank_id, unsigned int dpu_offset);

void dpu_try_get_results_and_log(unsigned int rank_id, unsigned int dpu_offset);

#endif /* __INTEGRATION_DPUS_H__ */
