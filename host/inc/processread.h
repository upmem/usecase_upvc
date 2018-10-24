#ifndef __PROCESSREAD_H__
#define __PROCESSREAD_H__

#include <stdint.h>
#include <stdio.h>
#include "genome.h"
#include "vartree.h"
#include "upvc.h"

// TODO: to be commented based on DL's feedbacks
// TODO: rename f1 f2
int process_read(genome_t *ref_genome,
                 int8_t *reads_buffer,
                 variant_tree_t **variant_list,
                 int *substitution_list,
                 int8_t *mapping_coverage,
                 FILE *f1,
                 FILE *f2,
                 int round,
                 times_ctx_t *times_ctx,
                 reads_info_t *reads_info);

#endif /* __PROCESSREAD_H__ */
