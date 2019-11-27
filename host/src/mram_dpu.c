/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mram_dpu.h"
#include "upvc.h"
#include "upvc_dpu.h"
#include "vmi.h"

#define MRAM_FILE_NAME_SIZE (24)

bool mram_copy_vmi(uint8_t *mram, vmi_t *vmi)
{
    vmi_read(vmi, mram, vmi->mem_size);
    return true;
}

static void make_mram_file_name(unsigned int dpu_id, char *name) { sprintf(name, "mram_%04u.bin", dpu_id); }

bool mram_save(uint8_t *mram, unsigned int dpu_id)
{
    char file_name[MRAM_FILE_NAME_SIZE];
    make_mram_file_name(dpu_id, file_name);
    FILE *f = fopen(file_name, "wb");
    assert(f != NULL);

    size_t write_size = MRAM_SIZE;
    size_t written = fwrite(mram, sizeof(uint8_t), write_size, f);
    if (written != write_size) {
        ERROR_EXIT(20, "BUG! fwrite in mram_save (%s) returned %lu instead of %lu - aborting!", file_name, written, write_size);
    }

    fclose(f);
    return true;
}

void mram_load(uint8_t *mram, unsigned int dpu_id)
{
    char file_name[MRAM_FILE_NAME_SIZE];
    make_mram_file_name(dpu_id, file_name);
    FILE *f = fopen(file_name, "rb");
    assert(f != NULL);

    __attribute__((unused)) size_t read_size = fread(mram, sizeof(uint8_t), MRAM_SIZE, f);

    fclose(f);
}

