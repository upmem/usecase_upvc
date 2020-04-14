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

#define MRAM_SIZE_AVAILABLE (MRAM_SIZE - MAX_DPU_REQUEST * sizeof(dpu_request_t) - MAX_DPU_RESULTS * sizeof(dpu_result_out_t))
static const int mram_size = MRAM_SIZE_AVAILABLE;
static vmi_t *vmis = NULL;
_Static_assert(MRAM_SIZE_AVAILABLE > 0, "Too many request and/or result compare to MRAM_SIZE");

static char *make_mram_file_name(unsigned int dpu_id)
{
    char *file_name;
    asprintf(&file_name, "mram_%04u.bin", dpu_id);
    assert(file_name != NULL);
    return file_name;
}

size_t mram_load(uint8_t **mram, unsigned int dpu_id)
{
    char * file_name = make_mram_file_name(dpu_id);
    FILE *f = fopen(file_name, "rb");
    assert(f != NULL);
    free(file_name);

    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    rewind(f);

    *mram = malloc(mram_size);
    assert(*mram != NULL);

    size_t read_size = fread(*mram, sizeof(uint8_t), size, f);
    assert(read_size == size);

    fclose(f);
    return read_size;
}

void init_vmis(unsigned int nb_dpu)
{
    vmis = (vmi_t *)calloc(nb_dpu, sizeof(vmi_t));
    assert(vmis != NULL);
    for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
        char *file_name = make_mram_file_name(dpuno);
        vmis[dpuno] = open(file_name, O_WRONLY | O_CREAT, S_IRUSR | S_IRGRP | S_IROTH);
        if (vmis[dpuno] == -1) {
            ERROR_EXIT(-1, "%s: open failed (%s)", file_name, strerror(errno));
        }
        free(file_name);
    }
}

void free_vmis(unsigned int nb_dpu)
{
    for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
        close(vmis[dpuno]);
    }
    free(vmis);
}

void write_vmi(unsigned int num_dpu, unsigned int num_ref, int8_t *nbr, dpu_result_coord_t coord)
{
    coords_and_nbr_t buffer;
    ssize_t buffer_size = sizeof(buffer);
    buffer.coord = coord;
    memcpy(&buffer.nbr, nbr, sizeof(buffer.nbr));

    vmi_t vmi = vmis[num_dpu];
    lseek(vmi, num_ref * sizeof(buffer), SEEK_SET);
    ssize_t write_size = write(vmi, &buffer, buffer_size);
    if (write_size != buffer_size) {
        ERROR_EXIT(-1, "%s: did not write the expected number of bytes (expected %li got %li, errno: %s)",
            make_mram_file_name(num_dpu), buffer_size, write_size, strerror(errno));
    }
}
