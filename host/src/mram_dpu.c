/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

char file_name[24] = { 0 };

char *make_mram_file_name(unsigned int dpu_id)
{
    sprintf(file_name, "mram_%04u.bin", dpu_id);
    return file_name;
}

size_t mram_load(uint8_t *mram, unsigned int dpu_id)
{
    FILE *f = fopen(make_mram_file_name(dpu_id), "rb");
    assert(f != NULL);

    size_t read_size = fread(mram, sizeof(uint8_t), MRAM_SIZE, f);

    fclose(f);
    return read_size;
}
