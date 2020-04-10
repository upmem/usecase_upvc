/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "common.h"
#include "mram_dpu.h"

char *make_mram_file_name(char *str, unsigned int dpu_id)
{
    sprintf(str, "mram_%04u.bin", dpu_id);
    return str;
}

size_t mram_load(uint8_t *mram, unsigned int dpu_id)
{
    char file_name[FILE_NAME_SIZE];
    FILE *f = fopen(make_mram_file_name(file_name, dpu_id), "rb");
    assert(f != NULL);

    size_t read_size = fread(mram, sizeof(uint8_t), MRAM_SIZE, f);

    fclose(f);
    return read_size;
}
