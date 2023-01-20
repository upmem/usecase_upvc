/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "debug.h"
#include "getread.h"
#include "upvc.h"

#define MAX_SEQ_SIZE (512)
#define MAX_BUF_SIZE (1024)

static int nb_reads[NB_READS_BUFFER];
static int8_t *reads_buffers[NB_READS_BUFFER];
static float *reads_quality_buffers[NB_READS_BUFFER];
static float quality_lookup_table[43];
static bool lookup = false;
#define PASS(pass_id) (pass_id % NB_READS_BUFFER)

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/**
 * @brief Parse the file "f" to get the next read in the file and its pair.
 *
 * @param f                     File to parse.
 * @param read1                 Output the next read in the file.
 * @param read2                 Output the pair of the next read in the file.
 * @param read_quality_factor   Output the quality factor of the read.
 *
 * @return The size of the read.
 */
static int get_seq_fast_AQ(FILE *f, int8_t *read1, int8_t *read2, float *read_quality_factor)
{
    static const int invnt[4] = { 2, 3, 0, 1 };
    int offset = 0;
    char comment[MAX_BUF_SIZE];
    char sequence_buffer[MAX_SEQ_SIZE];

    if (fgets(comment, MAX_BUF_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Sequence */
        return -1;
    }

    /* If the comment start with ">>14"
     * it means that we need the skip the first 14 characters of the read.
     */
    if (comment[1] == '>') {
          sscanf(&comment[2], "%d", &offset);
      }
      int i;
      for (i = 0; i < SIZE_READ - offset; i++) {
        read1[i] = (((int)sequence_buffer[i]) >> 1) & 3;
        read2[SIZE_READ - i - 1 - offset] = invnt[read1[i]];
    }
    for (; i < SIZE_READ; i++) {
        read1[i] = 0;
        read2[i] = 0;
    }

    if (comment[0] == '>') {
        return SIZE_READ;
    }
    if (fgets(comment, MAX_BUF_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Line with sequence quality information */
        return -1;
    }
    // store quality information
    for (i = 0; i < SIZE_READ - offset; i++) {
      int Q = sequence_buffer[i];
      read_quality_factor[i] = quality_lookup_table[Q-33];
    }
    for (; i < SIZE_READ; i++) {
      read_quality_factor[i] = 1.0f;
    }

    return SIZE_READ;
}

void get_reads(FILE *fpe1, FILE *fpe2, unsigned int pass_id)
{
    int nb_read = 0;
    pass_id = PASS(pass_id);

    if(!lookup) {

        for(uint8_t i = 33; i < 76; ++i) {
            quality_lookup_table[i - 33] = MAX(0.001f, 1.0f - pow(10.0f, -(i-33)/10.0f));
            /*printf("%d %f\n", i, quality_lookup_table[i - 33]);*/
        }
        lookup = true;
    }

    int8_t *reads_buffer = reads_buffers[pass_id];
    float* reads_quality_buffer = reads_quality_buffers[pass_id];
    if (reads_buffer == NULL) {
        LOG_INFO("allocating %lu for reads_buffer\n", MAX_READS_BUFFER*SIZE_READ*(1+sizeof(float)));
        reads_buffer = (int8_t *)malloc(MAX_READS_BUFFER * SIZE_READ);
        reads_quality_buffer = (float *)malloc(MAX_READS_BUFFER/2 * SIZE_READ * sizeof(float));
        assert(reads_buffer != NULL);
        assert(reads_quality_buffer != NULL);
        reads_buffers[pass_id] = reads_buffer;
        reads_quality_buffers[pass_id] = reads_quality_buffer;
    }

    while (nb_read < MAX_READS_BUFFER) {
        /*
        if (pass_id>=20) // FIXME : remove this condition used for debugging
            break;
            */
        if ((get_seq_fast_AQ(fpe1, &reads_buffer[(nb_read + 0) * SIZE_READ], &reads_buffer[(nb_read + 1) * SIZE_READ], 
                &reads_quality_buffer[nb_read/2 * SIZE_READ]) <= 0)
            || (get_seq_fast_AQ(fpe2, &reads_buffer[(nb_read + 2) * SIZE_READ], &reads_buffer[(nb_read + 3) * SIZE_READ], 
                &reads_quality_buffer[(nb_read/2 + 1) * SIZE_READ]) <= 0))
            break;
        nb_read += 4;
    }

    nb_reads[pass_id] = nb_read;
}

int get_reads_in_buffer(unsigned int pass_id) { return nb_reads[PASS(pass_id)]; }

int8_t *get_reads_buffer(unsigned int pass_id) { return reads_buffers[PASS(pass_id)]; }
float *get_reads_quality_buffer(unsigned int pass_id) { return reads_quality_buffers[PASS(pass_id)]; }

int get_input_info(FILE *f, size_t *read_size, size_t *nb_read)
{
    size_t size;
    char sequence_buffer[MAX_SEQ_SIZE];

    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Commentary */
        return -1;
    }
    if (fgets(sequence_buffer, MAX_SEQ_SIZE, f) == NULL) { /* Sequence */
        return -1;
    }
    size = strlen(sequence_buffer) - 1;

    fseek(f, 0, SEEK_END);
    size_t file_size = ftell(f);

    *read_size = size;
    *nb_read = file_size / (2 * size);

    rewind(f);
    return 0;
}
