/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genome.h"
#include "upvc.h"

#define MAX_BUF_SIZE (1024)

genome_t *get_genome(char *fasta_file, times_ctx_t *times_ctx)
{
    FILE *genome_file;
    size_t genome_file_size;
    double t1, t2;
    char genome_file_line[MAX_BUF_SIZE];
    genome_t *genome;

    printf("Read genome\n");
    t1 = my_clock();

    genome_file = fopen(fasta_file, "r");
    fseek(genome_file, 0, SEEK_END);
    genome_file_size = ftell(genome_file);
    rewind(genome_file);

    genome = (genome_t *)malloc(sizeof(genome_t));
    genome->data = (int8_t *)malloc(sizeof(int8_t) * genome_file_size);
    genome->fasta_file_name = strdup((const char *)fasta_file);
    genome->fasta_file_size = genome_file_size;
    genome->nb_seq = 0;

    {
        uint64_t current_data_idx = 0;
        while (fgets(genome_file_line, MAX_BUF_SIZE, genome_file) != NULL) {
            if (genome_file_line[0] == '>') { /* Commentary line with metadata */
                genome->pt_seq[genome->nb_seq] = current_data_idx;
                genome->len_seq[genome->nb_seq] = 0;
                genome_file_line[strlen(genome_file_line) - 1] = '\0';
                genome->seq_name[genome->nb_seq] = strndup(&genome_file_line[1], strlen(genome_file_line) - 1);
                genome->nb_seq++;
            } else {
                for (unsigned int i = 0; i < strlen(genome_file_line) - 1; i++) {
                    if (genome_file_line[i] != 'N') { /* A -> 0, C -> 1, G -> 3, T -> 2 */
                        genome->data[current_data_idx++] = (((int)genome_file_line[i]) >> 1) & 3;
                    } else {
                        genome->data[current_data_idx++] = 4;
                    }
                    genome->len_seq[genome->nb_seq - 1]++;
                }
            }
        }
    }
    fclose(genome_file);

    t2 = my_clock();
    times_ctx->get_genome = t2 - t1;

    printf(" - #seq: %d\n", genome->nb_seq);
    printf(" - time: %lf\n", times_ctx->get_genome);

    return genome;
}

void free_genome(genome_t *genome)
{
    for (uint32_t i = 0; i < genome->nb_seq; i++) {
        free(genome->seq_name[i]);
    }
    free(genome->fasta_file_name);
    free(genome->data);
    free(genome);
}
