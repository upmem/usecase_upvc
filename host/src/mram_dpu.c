/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _GNU_SOURCE
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"
#include "mram_dpu.h"
#include "upvc.h"

#define MRAM_FORMAT "mram_%04u.bin"
#define MRAM_SIZE_AVAILABLE (MRAM_SIZE - MAX_DPU_REQUEST * sizeof(dpu_request_t) - MAX_DPU_RESULTS * sizeof(dpu_result_out_t))
typedef struct {
    uint32_t size;
    uint8_t *buffer;
} vmi_t;
static const size_t mram_size = MRAM_SIZE_AVAILABLE;
static vmi_t *vmis = NULL;
_Static_assert(MRAM_SIZE_AVAILABLE > 0, "Too many request and/or result compare to MRAM_SIZE");

static char *make_mram_file_name(unsigned int dpu_id)
{
    char *file_name;
    char *index_folder = get_index_folder();
    assert(asprintf(&file_name, "%s" MRAM_FORMAT, index_folder, dpu_id) > 0);
    assert(file_name != NULL);
    free(index_folder);
    return file_name;
}

size_t mram_load(uint8_t **mram, unsigned int dpu_id)
{
    char *file_name = make_mram_file_name(dpu_id);
    FILE *f = fopen(file_name, "rb");
    CHECK_FILE(f, file_name);
    free(file_name);

    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    rewind(f);

    assert(mram_size >= size);
    *mram = malloc(mram_size);
    assert(*mram != NULL);

    size_t read_size = fread(*mram, sizeof(uint8_t), size, f);
    assert(read_size == size);

    fclose(f);
    return read_size;
}

void init_vmis(unsigned int nb_dpu, distribute_index_t *table)
{
    vmis = (vmi_t *)calloc(nb_dpu, sizeof(vmi_t));
    assert(vmis != NULL);

    for (unsigned int i = 0; i < nb_dpu; i++) {
        vmis[i].size = table[i].size * sizeof(coords_and_nbr_t);
        assert(vmis[i].size < MRAM_SIZE);
        vmis[i].buffer = (uint8_t *)malloc(vmis[i].size);
        assert(vmis[i].buffer != NULL);
    }
}

void free_vmis(unsigned int nb_dpu)
{
    for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
        char *file_name = make_mram_file_name(dpuno);
        FILE *f = fopen(file_name, "w");
        CHECK_FILE(f, file_name);
        free(file_name);
        xfer_file(vmis[dpuno].buffer, vmis[dpuno].size, f, xfer_write);
        fclose(f);
        free(vmis[dpuno].buffer);
    }

    free(vmis);
}

void write_vmi(unsigned int num_dpu, unsigned int num_ref, int8_t *nbr, dpu_result_coord_t coord)
{
    uint32_t offset = sizeof(coords_and_nbr_t) * num_ref;
    assert(offset < vmis[num_dpu].size);
    coords_and_nbr_t *buffer = (coords_and_nbr_t *)&vmis[num_dpu].buffer[offset];
    buffer->coord = coord;
    memcpy(&buffer->nbr, nbr, sizeof(buffer->nbr));
}
