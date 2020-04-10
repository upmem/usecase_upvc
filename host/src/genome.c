/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "genome.h"
#include "upvc.h"

#define MAX_BUF_SIZE (1024)

static genome_t genome;

genome_t *genome_get() { return &genome; }

void genome_init(char *fasta_file)
{
    char genome_file_line[MAX_BUF_SIZE];

    double start_time = my_clock();

    printf("%s:\n", __func__);

    FILE *genome_file = fopen(fasta_file, "r");
    assert(genome_file != NULL);

    fseek(genome_file, 0, SEEK_END);
    genome.fasta_file_size = ftell(genome_file);
    rewind(genome_file);

    genome.data = (int8_t *)malloc(sizeof(int8_t) * genome.fasta_file_size);
    genome.mapping_coverage = (int32_t *)calloc(sizeof(int32_t), genome.fasta_file_size);
    assert(genome.data != NULL && genome.mapping_coverage != NULL);
    genome.nb_seq = 0;

    uint64_t current_data_idx = 0;
    while (fgets(genome_file_line, MAX_BUF_SIZE, genome_file) != NULL) {
        if (genome_file_line[0] == '>') { /* Commentary line with metadata */
            genome.pt_seq[genome.nb_seq] = current_data_idx;
            genome.len_seq[genome.nb_seq] = 0;
            genome_file_line[strlen(genome_file_line) - 1] = '\0';
            genome.seq_name[genome.nb_seq] = strndup(&genome_file_line[1], strlen(genome_file_line) - 1);
            genome.nb_seq++;
        } else {
            for (unsigned int i = 0; i < strlen(genome_file_line) - 1; i++) {
                if (genome_file_line[i] != 'N') { /* A -> 0, C -> 1, G -> 3, T -> 2 */
                    genome.data[current_data_idx++] = (((int)genome_file_line[i]) >> 1) & 3;
                } else {
                    genome.data[current_data_idx++] = 4;
                }
                genome.len_seq[genome.nb_seq - 1]++;
            }
        }
    }
    fclose(genome_file);

    printf("\t#seq: %d\n", genome.nb_seq);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void genome_free()
{
    for (uint32_t i = 0; i < genome.nb_seq; i++) {
        free(genome.seq_name[i]);
    }
    free(genome.data);
    free(genome.mapping_coverage);
}
