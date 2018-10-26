/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_MDPU_H__
#define __INTEGRATION_MDPU_H__

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
 * 64K results of 8 bytes.
 */
#define MAX_DPU_RESULTS (1 << 16)
#define RESULT_AREA_LEN (MAX_DPU_RESULTS * 16)

/**
 * @brief Total size of the input area, in bytes: the whole MRAM minus the result area MINUS the swap area used
 * by the DPU internally.
 */
#define MRAM_INPUT_SIZE (MRAM_SIZE - RESULT_AREA_LEN - (16 * 1024 * 16))

/**
 * @brief A snapshot of the MRAM.
 *
 * @var io_offs   Offset to the I/O (reads and results) area in MRAM.
 * @var usage     Track the memory usage.
 * @var nb_nbr    Total number of neighbours stored in this MRAM.
 * @var nbr_offs  Offsets to the persistent neighbourhood data.
 * @var nbr_len   Length of a neighbour, in bytes.
 * @var magic     Magic number, for debugging purpose.
 */
typedef struct {
    uint32_t io_offs;
    uint32_t usage;
    uint32_t nb_nbr;
    uint32_t nbr_offs;
    uint32_t nbr_len;
    uint32_t magic;
} mram_info_t;

/**
 * @brief Creates a new MRAM image.
 *
 * The returned value represents one chunk of 64MBs to store an MRAM mapping.
 *
 * @param reads_info  Information on the size of the seed and the neighbour.
 *
 * @return The created image
 */
mram_info_t *mram_create(reads_info_t *reads_info);

/**
 * @brief Gets rid of an MRAM image.
 *
 * @param mram  The MRAM image, not usable anymore.
 */
void mram_free(mram_info_t *mram);

/**
 * @brief Initializes a mapping.
 *
 * @param mram        The target MRAM image, all fields will be reset.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
void mram_reset(mram_info_t *mram, reads_info_t *reads_info);

/**
 * @param mram  An MRAM image.
 *
 * @return The neighbour area, including coordinates.
 */
static inline long *mram_neighbours_area(mram_info_t *mram) {
    /* Suspicious pointer cast, but okay because nbr_offs is aligned on longs. */
    return (long *) (((uint8_t *) mram) + (mram->nbr_offs));
}

/**
 * @brief Puts a VMI file directly into the neighbour and coordinates area of an MRAM file.
 *
 * Creates the MRAM header accordingly.
 *
 * @param mram        The target MRAM image.
 * @param vmi         The VMI.
 * @param nb_nbr      How many neighbours are registered.
 * @param reads_info  Information on the size of the seed and the neighbour..
 *
 * @return Whether the VMI is copied (is false if the MRAM is full).
 */
bool mram_copy_vmi(mram_info_t *mram, vmi_t *vmi, unsigned int nb_nbr, reads_info_t *reads_info);

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

#endif /* __INTEGRATION_MDPU_H__ */
