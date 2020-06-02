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
#include "parse_args.h"
#include "upvc.h"

#include <dpu.h>

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

static struct dpu_set_t dpu_set;
static struct dpu_set_t *dpus;
static uint32_t nb_dpu_set;
static struct dpu_symbol_t mram_symbol = { .address = 0x08000000, .size = MRAM_SIZE };

void init_vmis(unsigned int nb_dpu, distribute_index_t *table)
{
    check_ulimit_n(nb_dpu * 2);

    vmis = (vmi_t *)calloc(nb_dpu, sizeof(vmi_t));
    assert(vmis != NULL);

    dpu_error_t err = !DPU_OK;
    if (get_index_with_dpus()) {
        err = dpu_alloc(DPU_ALLOCATE_ALL, "backend=hw", &dpu_set);
    }
    if (err == DPU_OK) {
        struct dpu_set_t dpu;
        uint32_t each_dpu;
        DPU_ASSERT(dpu_get_nr_dpus(dpu_set, &nb_dpu_set));
        printf("\t\tUsing %u DPUs to help indexing\n", nb_dpu_set);
        if (nb_dpu_set > nb_dpu) {
            nb_dpu_set = nb_dpu;
        }
        dpus = (struct dpu_set_t *)malloc(sizeof(struct dpu_set_t) * nb_dpu_set);
        assert(dpus != NULL);
        DPU_FOREACH (dpu_set, dpu, each_dpu) {
            if (each_dpu >= nb_dpu_set) {
                break;
            }
            dpus[each_dpu] = dpu;
        }
    } else {
        nb_dpu_set = 0;
    }

    for (unsigned int i = 0; i < nb_dpu; i++) {
        vmis[i].size = table[i].size * sizeof(coords_and_nbr_t);
        assert(vmis[i].size < MRAM_SIZE);
        if (i >= nb_dpu_set) {
            vmis[i].buffer = (uint8_t *)malloc(vmis[i].size);
            assert(vmis[i].buffer != NULL);
        }
    }
}

void free_vmis(unsigned int nb_dpu)
{
    uint8_t *tmp_mram = (uint8_t *)malloc(MRAM_SIZE);
    assert(tmp_mram != NULL);
    for (unsigned int dpuno = 0; dpuno < nb_dpu_set; dpuno++) {
        char *file_name = make_mram_file_name(dpuno);
        FILE *f = fopen(file_name, "w");
        CHECK_FILE(f, file_name);
        free(file_name);
        DPU_ASSERT(dpu_copy_from_symbol(dpus[dpuno], mram_symbol, 0, tmp_mram, vmis[dpuno].size));
        xfer_file(tmp_mram, vmis[dpuno].size, f, xfer_write);
        fclose(f);
    }
    for (unsigned int dpuno = nb_dpu_set; dpuno < nb_dpu; dpuno++) {
        char *file_name = make_mram_file_name(dpuno);
        FILE *f = fopen(file_name, "w");
        CHECK_FILE(f, file_name);
        free(file_name);
        xfer_file(vmis[dpuno].buffer, vmis[dpuno].size, f, xfer_write);
        fclose(f);
        free(vmis[dpuno].buffer);
    }

    free(tmp_mram);
    free(vmis);
}

void write_vmi(unsigned int num_dpu, unsigned int num_ref, coords_and_nbr_t *coords_and_nbr)
{
    uint32_t offset = sizeof(coords_and_nbr_t) * num_ref;
    assert(offset < vmis[num_dpu].size);
    if (num_dpu < nb_dpu_set) {
        DPU_ASSERT(dpu_copy_to_symbol(dpus[num_dpu], mram_symbol, offset, coords_and_nbr, sizeof(*coords_and_nbr)));
        return;
    }
    coords_and_nbr_t *buffer = (coords_and_nbr_t *)&vmis[num_dpu].buffer[offset];
    memcpy(buffer, coords_and_nbr, sizeof(*coords_and_nbr));
}
