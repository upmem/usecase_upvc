/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upvc.h"
#include "common.h"

#define MAX_SEQ_SIZE (512)
#define MAX_BUF_SIZE (1024)

static int nb_reads[NB_READS_BUFFER];
static int8_t *reads_buffers[NB_READS_BUFFER];

/**
 * @brief Parse the file "f" to get the next read in the file and its pair.
 *
 * @param f          File to parse.
 * @param read1      Output the next read in the file.
 * @param read2      Output the pair of the next read in the file.
 *
 * @return The size of the read.
 */
static int get_seq_fast_AQ(FILE *f, int8_t *read1, int8_t *read2)
{
    int offset = 0;
    int size_read = SIZE_READ;
    char comment[MAX_BUF_SIZE];
    char sequence_buffer[MAX_SEQ_SIZE];
    int invnt[4] = { 2, 3, 0, 1 };

    if (fgets(comment, MAX_BUF_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Sequence */
        return -1;
    }

    /* If the comment start with ">>12"
     * it means that we need the skip the first 12 characters of the read.
     */
    if (comment[1] == '>') {
        sscanf(&comment[2], "%d", &offset);
    }
    for (int i = 0; i < size_read - offset; i++) {
        read1[i] = (((int)sequence_buffer[i]) >> 1) & 3;
        read2[size_read - i - 1 - offset] = invnt[read1[i]];
    }
    for (int i = size_read - offset; i < size_read; i++) {
        read1[i] = 0;
        read2[i] = 0;
    }

    if (comment[0] == '>') {
        return size_read;
    }
    if (fgets(comment, MAX_BUF_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Line with sequence quality information (unused) */
        return -1;
    }
    return size_read;
}

void get_reads(FILE *fpe1, FILE *fpe2, unsigned int pass_id)
{
    int nb_read = 0;
    int size_read = SIZE_READ;
    pass_id %= NB_READS_BUFFER;

    int8_t *reads_buffer = reads_buffers[pass_id];
    if (reads_buffer == NULL) {
        reads_buffer = (int8_t *)malloc(MAX_READS_BUFFER * SIZE_READ);
        assert(reads_buffer != NULL);
        reads_buffers[pass_id] = reads_buffer;
    }

    while (nb_read < MAX_READS_BUFFER) {
        if ((get_seq_fast_AQ(fpe1, &reads_buffer[(nb_read + 0) * size_read], &reads_buffer[(nb_read + 1) * size_read]) <= 0)
            || (get_seq_fast_AQ(fpe2, &reads_buffer[(nb_read + 2) * size_read], &reads_buffer[(nb_read + 3) * size_read]) <= 0))
            break;
        nb_read += 4;
    }

    nb_reads[pass_id] = nb_read;
}

int get_reads_in_buffer(unsigned int pass_id) { return nb_reads[pass_id % NB_READS_BUFFER]; }

int8_t *get_reads_buffer(unsigned int pass_id) { return reads_buffers[pass_id % NB_READS_BUFFER]; }

int get_read_size(char *input_pe1_file)
{
    FILE *f;
    size_t size;
    char sequence_buffer[MAX_SEQ_SIZE];

    f = fopen(input_pe1_file, "r");
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Sequence */
        return -1;
    }
    size = strlen(sequence_buffer) - 1;
    size &= ~((int)0x3);

    fclose(f);
    return (int)size;
}
