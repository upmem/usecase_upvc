#ifndef __VCF_H__
#define __VCF_H__

#include <stdint.h>
#include "genome.h"
#include "vartree.h"
#include "upvc.h"

/*
 * Produce the final VCF file
 */
int create_vcf(char *chromosome_name,
                genome_t *ref_genome,
                variant_tree_t **variant_list,
                int *substitution_list,
                int8_t *mapping_coverage,
                times_ctx_t *times_ctx);

#endif /*__VCF_H__ */
