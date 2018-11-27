/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdint.h>

#define MAX_SEQ_GEN (1000)

/**
 * @brief Reference genome information.
 *
 * @var data             Byte list of the genome.
 * @var nb_seq           Numbers of sequences composing the genome.
 * @var pt_seq           Index in data of the start of each sequence.
 * @var len_seq          Sequence lengths.
 * @var fasta_file_size  Size of the original file.
 * @var seq_name         Sequence names.
 * @var fasta_file_name  Name of the original file.
 */
typedef struct {
        int8_t *data;
        int nb_seq;
        int pt_seq[MAX_SEQ_GEN];
        int len_seq[MAX_SEQ_GEN];
        long fasta_file_size;
        char *seq_name[MAX_SEQ_GEN];
        char *fasta_file_name;
} genome_t;

#include "upvc.h"

/**
 * @brief Function that create the structure holding the information on the reference genome.
 *
 * @param filename_prefix  File to be parsed to find the genome information.
 * @param times_ctx        Times information for the whole application.
 *
 * @return The structure holding the informatino on the reference genome.
 */
genome_t *get_genome(char* filename_prefix, times_ctx_t *times_ctx);

/**
 * @brief Free a "genome_t" structure.
 *
 * @param genome  The "genome_t" structure to be freed.
 */
void free_genome(genome_t *genome);

#endif /* __GENOME_H__ */
