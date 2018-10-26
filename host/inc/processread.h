/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __PROCESSREAD_H__
#define __PROCESSREAD_H__

#include <stdint.h>
#include <stdio.h>
#include "genome.h"
#include "vartree.h"
#include "upvc.h"

/**
 * @brief Add each read that matched into the variant list.
 * Reads that didn't match are to be print into file to be reprocess.
 *
 * @param ref_genome         Reference genome.
 * @param reads_buffer       Input buffer of reads.
 * @param variant_list       Output list of the variant found in process_read.
 * @param substitution_list  List of substitution.
 * @param mapping_coverage   Number of match for every position of the reference genome.
 * @param fpe1               First file of Pair-end read containing the read that didn't match.
 * @param fpe2               Second file of Pair-end read containing the read that didn't match.
 * @param round              Round number of the process_read.
 * @param nb_dpu             Number of DPUs used to compute.
 * @param times_ctx          Times information for the whole application.
 * @param reads_info         Information on the size of the seed and the neighbour.
 *
 * @return The number of reads that matched.
 */
int process_read(genome_t *ref_genome,
                 int8_t *reads_buffer,
                 variant_tree_t **variant_list,
                 int *substitution_list,
                 int8_t *mapping_coverage,
                 FILE *fpe1,
                 FILE *fpe2,
                 int round,
                 int nb_dpu,
                 times_ctx_t *times_ctx,
                 reads_info_t *reads_info);

#endif /* __PROCESSREAD_H__ */
