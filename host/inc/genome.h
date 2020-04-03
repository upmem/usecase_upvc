/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdint.h>

#define MAX_SEQ_GEN (24) // number of chromosome

/**
 * @brief Reference genome information.
 *
 * @var data             Byte list of the genome.
 * @var nb_seq           Numbers of sequences composing the genome.
 * @var pt_seq           Index in data of the start of each sequence.
 * @var len_seq          Sequence lengths.
 * @var fasta_file_size  Size of the original file.
 * @var seq_name         Sequence names.
 */
typedef struct {
    int8_t *data;
    int32_t *mapping_coverage;
    int *substitution_list;
    uint32_t nb_seq;
    uint64_t pt_seq[MAX_SEQ_GEN];
    uint64_t len_seq[MAX_SEQ_GEN];
    uint64_t fasta_file_size;
    char *seq_name[MAX_SEQ_GEN];
} genome_t;

void genome_init(char *filename_prefix);

void genome_free();

genome_t *genome_get();

#endif /* __GENOME_H__ */
