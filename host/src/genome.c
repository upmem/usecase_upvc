/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "index.h"
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "genome.h"
#include "parse_args.h"
#include "upvc.h"

#define MAX_BUF_SIZE (1024)
#define GENOME_BINARY "genome.bin"
#define GENOME_VERSION (1)
#define GENOME_MAGIC (0x9e503e)

static genome_t genome = { .magic = GENOME_MAGIC, .version = GENOME_VERSION };

genome_t *genome_get() { return &genome; }

void genome_load()
{
    double start_time = my_clock();

    printf("%s:\n", __func__);
    char filename[FILENAME_MAX];
    char *index_folder = get_index_folder();
    sprintf(filename, "%s" GENOME_BINARY, index_folder);
    free(index_folder);
    FILE *f = fopen(filename, "r");
    CHECK_FILE(f, filename);

    xfer_file((uint8_t *)&genome, sizeof(genome_t), f, xfer_read);

    assert(genome.magic == GENOME_MAGIC
        && "Wrong header, make sure you have generated your MRAMs with the same version of UPVC that you are using.");
    assert(genome.version == GENOME_VERSION && "Could not load a genome generated with a different version of UPVC.");

    genome.data = (int8_t *)malloc(sizeof(int8_t) * genome.fasta_file_size);
    genome.mapping_coverage = (int32_t *)calloc(sizeof(int32_t), genome.fasta_file_size);
    assert(genome.data != NULL && genome.mapping_coverage != NULL);

    xfer_file((uint8_t *)genome.data, genome.fasta_file_size, f, xfer_read);

    fclose(f);

    printf("\t#seq: %d\n", genome.nb_seq);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void genome_create()
{
    char *prefix = get_input_path();
    char genome_file_line[MAX_BUF_SIZE];

    double start_time = my_clock();

    printf("%s:\n", __func__);

    char filename[FILENAME_MAX];
    sprintf(filename, "%s.fasta", prefix);
    FILE *genome_file = fopen(filename, "r");
    CHECK_FILE(genome_file, filename);

    fseek(genome_file, 0, SEEK_END);
    genome.fasta_file_size = ftell(genome_file);
    rewind(genome_file);

    genome.data = (int8_t *)malloc(sizeof(int8_t) * genome.fasta_file_size);
    assert(genome.data != NULL);
    genome.nb_seq = 0;
    genome.mapping_coverage = NULL;

    uint64_t current_data_idx = 0;
    while (fgets(genome_file_line, MAX_BUF_SIZE, genome_file) != NULL) {
        if (genome_file_line[0] == '>') { /* Commentary line with metadata */
            genome.pt_seq[genome.nb_seq] = current_data_idx;
            genome.len_seq[genome.nb_seq] = 0;
            genome_file_line[strlen(genome_file_line) - 1] = '\0';
            memcpy(genome.seq_name[genome.nb_seq], &genome_file_line[1],
                strlen(&genome_file_line[1]) > MAX_SEQ_NAME_SIZE ? MAX_SEQ_NAME_SIZE : strlen(&genome_file_line[1]));
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

    char *index_folder = get_index_folder();
    sprintf(filename, "%s" GENOME_BINARY, index_folder);
    free(index_folder);
    FILE *genome_binary = fopen(filename, "w");
    CHECK_FILE(genome_binary, filename);
    xfer_file((uint8_t *)&genome, sizeof(genome_t), genome_binary, xfer_write);
    xfer_file((uint8_t *)genome.data, genome.fasta_file_size, genome_binary, xfer_write);
    fclose(genome_binary);

    printf("\t#seq: %d\n", genome.nb_seq);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void genome_free()
{
    free(genome.data);
    free(genome.mapping_coverage);
}
