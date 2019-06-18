/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
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

static void mram_reset(mram_info_t *mram, reads_info_t *reads_info)
{
    mram->total_nbr_size = 0;
    mram->nb_nbr = 0;
    mram->nbr_len = (uint32_t)reads_info->size_neighbour_in_bytes;
    mram->delta = reads_info->delta_neighbour_in_bytes;
}

bool mram_copy_vmi(mram_info_t *mram, vmi_t *vmi, unsigned int nb_nbr, reads_info_t *reads_info)
{
    if (vmi->mem_size + sizeof(mram_info_t) >= MRAM_SIZE) {
        ERROR_EXIT(18, "MRAM size exceeded when copying %u bytes of neighbourhood", (unsigned int)vmi->mem_size);
    }

    mram_reset(mram, reads_info);
    mram->nb_nbr = nb_nbr;

    vmi_read(vmi, ((uint8_t *)mram) + ALIGN_DPU(sizeof(mram_info_t)), vmi->mem_size);
    mram->total_nbr_size = ALIGN_DPU(vmi->mem_size);
    return true;
}

static void make_mram_file_name(unsigned int dpu_id, char *name) { sprintf(name, "mram_%04u.bin", dpu_id); }

bool mram_save(mram_info_t *mram, unsigned int dpu_id)
{
    char file_name[MRAM_FILE_NAME_SIZE];
    make_mram_file_name(dpu_id, file_name);
    FILE *f = fopen(file_name, "wb");
    assert(f != NULL);

    size_t write_size = ALIGN_DPU(mram->total_nbr_size + sizeof(mram_info_t));
    size_t written = fwrite(mram, sizeof(uint8_t), write_size, f);
    if (written != write_size) {
        ERROR_EXIT(20, "BUG! fwrite in mram_save (%s) returned %lu instead of %lu - aborting!", file_name, written, write_size);
    }

    fclose(f);
    return true;
}

static void mram_load_size(mram_info_t *mram, unsigned int dpu_id, unsigned int size)
{
    char file_name[MRAM_FILE_NAME_SIZE];
    make_mram_file_name(dpu_id, file_name);
    FILE *f = fopen(file_name, "rb");
    assert(f != NULL);

    __attribute__((unused)) size_t read_size = fread(mram, sizeof(uint8_t), size, f);

    fclose(f);
}
void mram_load(mram_info_t *mram, unsigned int dpu_id) { mram_load_size(mram, dpu_id, MRAM_SIZE); }

void mram_load_info(mram_info_t *mram, unsigned int dpu_id) { mram_load_size(mram, dpu_id, sizeof(mram_info_t)); }
