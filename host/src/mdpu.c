/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <string.h>
#include <stdlib.h>

#include "mdpu.h"
#include "vmi.h"
#include "mdpu.h"
#include "upvc.h"

#define MRAM_FILE_NAME_SIZE (24)

mram_info_t *mram_create(reads_info_t *reads_info)
{
        mram_info_t *result = (mram_info_t *) malloc(MRAM_SIZE);
        if (result == NULL) {
                ERROR_EXIT(17, "*** malloc failed to create a new MRAM image!");
        }
        mram_reset(result, reads_info);
        return result;
}

void mram_free(mram_info_t *mram)
{
        free(mram);
}

void mram_reset(mram_info_t *mram, reads_info_t *reads_info)
{
        mram->io_offs = 0;
        mram->magic = 0x42424242;
        mram->nbr_offs = 0;
        mram->nb_nbr = 0;
        mram->nbr_len = (uint32_t) reads_info->size_neighbour_in_bytes;
        mram->usage = 6 * sizeof(uint32_t);
}

bool mram_copy_vmi(mram_info_t *mram, vmi_t *vmi, unsigned int nb_nbr, reads_info_t *reads_info)
{
        if (vmi->mem_size + sizeof(mram_info_t) >= MRAM_SIZE) {
                ERROR_EXIT(18, "MRAM size exceeded when copying %u bytes of neighbourhood", (unsigned int) vmi->mem_size);
        }

        mram_reset(mram, reads_info);
        mram->nbr_offs = (mram->usage + 7) & ~0x7;
        mram->nb_nbr = nb_nbr;

        vmi_read(vmi, 0, ((uint8_t *) mram) + sizeof(mram_info_t), vmi->mem_size);
        mram->io_offs = mram->usage = (sizeof(mram_info_t) + vmi->mem_size + 7) & ~7;
        return true;
}


static void make_mram_file_name(unsigned int dpu_id, char *name)
{
        sprintf(name, "mram_%04u.bin", dpu_id);
}

bool mram_save(mram_info_t *mram, unsigned int dpu_id) {
        char file_name[MRAM_FILE_NAME_SIZE];
        make_mram_file_name(dpu_id, file_name);
        FILE *f = fopen(file_name, "wb");
        if (f == NULL) {
                ERROR_EXIT(19, "FATAL: could not create file '%s' for writing - aborting!", file_name);
        }

        size_t written = fwrite(mram, sizeof(uint8_t), mram->usage, f);
        if (written != mram->usage) {
                ERROR_EXIT(20, "BUG! fwrite in mram_save (%s) returned %lu instead of %u - aborting!", file_name, written,
                           mram->usage);
        }

        fclose(f);
        return true;
}

void mram_load(mram_info_t *mram, unsigned int dpu_id) {
        char file_name[MRAM_FILE_NAME_SIZE];
        make_mram_file_name(dpu_id, file_name);
        FILE *f = fopen(file_name, "rb");
        if (f == NULL) {
                ERROR_EXIT(21, "could not load MRAM file '%s'", file_name);
        }

        fread(mram, sizeof(uint8_t), MRAM_SIZE, f);
}
