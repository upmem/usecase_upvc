/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __VCF_H__
#define __VCF_H__

#include <stdint.h>
#include "genome.h"
#include "vartree.h"
#include "upvc.h"

/**
 * @brief Produce the final VCF file.
 *
 * @param chromosome_name    Name of the chromosone.
 * @param ref_genome         Reference genome.
 * @param variant_list       List (tree) of result variant.
 * @param substitution_list  List of substitution.
 * @param mapping_coverage   Number of match for every position of the reference genome.
 * @param times_ctx          Times information for the whole application.
 */
void create_vcf(char *chromosome_name,
                genome_t *ref_genome,
                variant_tree_t **variant_list,
                int *substitution_list,
                int8_t *mapping_coverage,
                times_ctx_t *times_ctx);

#endif /*__VCF_H__ */
