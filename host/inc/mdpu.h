/*
 * Copyright (c) 2014-2018 - uPmem
 */
#ifndef INTEGRATION_MDPU_H
#define INTEGRATION_MDPU_H

/**
 * @brief Mapping of the DPU MRAMs.
 *
 * Defines the structures representing the DPU MRAMs on both the host and DPU side.
 */
#include <stdint.h>
#include <stdbool.h>

#include "dpus.h"
#include "vmi.h"
#include "upvc.h"

#define MRAM_SIZE (64 << 20)

/**
 * @brief The MRAM is divided in an "input" area and a small area to store DPU results: the output area stores at most
 * 64K results of 8 bytes
 */
#define MAX_DPU_RESULTS (1 << 16)
#define RESULT_AREA_LEN (MAX_DPU_RESULTS * 16)

/**
 * @brief Total size of the input area, in bytes: the whole MRAM minus the result area MINUS the swap area used
 * by the DPU internally.
 */
#define MRAM_INPUT_SIZE (MRAM_SIZE - RESULT_AREA_LEN - (16 * 1024 * 16))

/**
 * @brief A snapshot of the MRAM
 * @var io_offs      offset to the I/O (reads and results) area in MRAM
 * @var usage        track the memory usage
 * @var nr_nbr       total number of neighbors stored in this MRAM
 * @var nbr_offs     offsets to the persistent neighborhood data
 * @var nbr_len      length of a neighbor, in bytes
 * @var magic        magic number, for debugging purpose
 */
typedef struct {
    uint32_t io_offs;
    uint32_t usage;
    uint32_t nr_nbr;
    uint32_t nbr_offs;
    uint32_t nbr_len;
    uint32_t magic;
} mram_info_t;

/**
 * @brief creates a new MRAM image
 *
 * The returned value represents one chunk of 64MBs to store an MRAM mapping.
 *
 * @param reads_info are the informations about the seed and neighbour sizes
 * @return the created image
 */
mram_info_t *mram_create(reads_info_t *reads_info);

/**
 * @brief gets rid of an MRAM image
 * @param mram the MRAM image, not usable anymore
 */
void mram_free(mram_info_t *mram);

/**
 * @brief initializes a mapping
 * @param mram the target MRAM image, all fields will be reset
 * @param reads_info are the informations about the seed and neighbour sizes
 */
void mram_reset(mram_info_t *mram, reads_info_t *reads_info);

/**
 * @param mram an MRAM image
 * @return the neighbor area, including coordinates
 */
static inline long *mram_neighbors_area(mram_info_t *mram) {
    // Suspicious pointer cast, but okay because nbr_offs is aligned on longs.
    return (long *) (((uint8_t *) mram) + (mram->nbr_offs));
}

/**
 * @brief puts a VMI file directly into the neighbor and coordinates area of an MRAM file.
 *
 * Creates the MRAM header accordingly.
 *
 * @param mram the target MRAM image
 * @param vmi  the VMI
 * @param nr_nbr how many neighbors are registered
 * @param reads_info are the informations about the seed and neighbour sizes
 * @return whether the VMI is copied (is false if the MRAM is full)
 */
bool mram_copy_vmi(mram_info_t *mram, vmi_t *vmi, unsigned int nr_nbr, reads_info_t *reads_info);

/**
 * @brief makes a snapshot of an MRAM image into a file.
 *
 * The output file is called "mram_XXXX.bin", where XXX is the DPU number.
 *
 * @param mram   the MRAM image
 * @param dpu_id index of the DPU
 * @return false if the snapshot could not be created
 */
bool mram_save(mram_info_t *mram, unsigned int dpu_id);

/**
 * @brief Loads an MRAM file into an MRAM image.
 * @param mram the MRAM image
 * @param dpu_id index of the DPU
 */
void mram_load(mram_info_t *mram, unsigned int dpu_id);

#endif //INTEGRATION_MDPU_H
