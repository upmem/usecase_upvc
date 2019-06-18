/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __PROCESSREAD_H__
#define __PROCESSREAD_H__

#include "common.h"
#include "genome.h"
#include "upvc.h"
#include "vartree.h"
#include <stdint.h>
#include <stdio.h>

/**
 * @brief Add each read that matched into the variant list.
 * Reads that didn't match are to be print into file to be reprocess.
 *
 * @param ref_genome         Reference genome.
 * @param reads_buffer       Input buffer of reads.
 * @param variant_list       Output list of the variant found in process_read.
 * @param substitution_list  List of substitution.
 * @param mapping_coverage   Number of match for every position of the reference genome.
 * @param result_tab         Buffer to store and sort results from every DPUs.
 * @param result_tab_nb_read Number of valid read in each index of the result tab.
 * @param fpe1               First file of Pair-end read containing the read that didn't match.
 * @param fpe2               Second file of Pair-end read containing the read that didn't match.
 * @param round              Round number of the process_read.
 * @param round              Current round to process.
 * @param times_ctx          Times information for the whole application.
 * @param reads_info         Information on the size of the seed and the neighbour.
 *
 * @return The number of reads that matched.
 */
unsigned int process_read(genome_t *ref_genome, int8_t *reads_buffer, variant_tree_t **variant_list, int *substitution_list,
    int8_t *mapping_coverage, dpu_result_out_t *result_tab, unsigned int result_tab_nb_read, FILE *fpe1, FILE *fpe2, int round,
    times_ctx_t *times_ctx, reads_info_t *reads_info);

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

#endif /* __PROCESSREAD_H__ */
