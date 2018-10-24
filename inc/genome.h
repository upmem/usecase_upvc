#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdint.h>
#include "upvc.h"

#define MAX_SEQ_GEN (1000)

typedef struct {
        int8_t *data;                   /* Byte list of the genome                        */
        int nb_seq;                     /* Numbers of sequences composing the genome      */
        int pt_seq[MAX_SEQ_GEN];        /* Index in data of the start of each sequence    */
        int len_seq[MAX_SEQ_GEN];       /* Sequence lengths                               */
        long fasta_file_size;           /* Size of the original file                      */
        char *seq_name[MAX_SEQ_GEN];    /* Sequence names                                 */
        char *fasta_file_name;          /* Name of the original file                      */
} genome_t;

/*
 * Parse the file starting with "filename_prefix".
 * Create the "genome_t" structure from it and return it.
 */
genome_t *get_genome(char* filename_prefix, times_ctx_t *times_ctx);
/*
 * Free a "genome_t" structure
 */
void free_genome(genome_t *G);

#endif /* __GENOME_H__ */
