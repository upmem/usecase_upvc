/**
 * Copyright 2020 - Dominique Lavenier & UPMEM
 */

#ifndef __ACCUMULATEREAD_H__
#define __ACCUMULATEREAD_H__

#include "common.h"
#include "upvc.h"

/**
 * @brief Accumulate read in result_tab in ordre to be process when the pass will have been execute on each DPUs.
 *
 * @param result_tab          Buffer to store and sort results from every DPUs of a pass.
 * @param result_tab_nb_read  Number of valid read in each index of the result tab.
 * @param dpu_offset          Offset in the virtual DPUs of the first DPU of this run.
 * @param times_ctx           Times information for the whole application.
 *
 * @return The number of reads accumulated.
 */
unsigned int accumulate_read(
    dpu_result_out_t **result_tab, unsigned int *result_tab_nb_read, unsigned int dpu_offset, times_ctx_t *times_ctx);

#endif /* __ACCUMULATEREAD_H__ */
