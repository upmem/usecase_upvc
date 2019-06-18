/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_MDPU_H__
#define __INTEGRATION_MDPU_H__

/**
 * @brief Mapping of the DPU MRAMs.
 *
 * Defines the structures representing the DPU MRAMs on both the host and DPU side.
 */
#include <stdbool.h>

#include "index.h"
#include "upvc.h"
#include "vmi.h"

#include "common.h"

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
 *
 * @return false if the snapshot could not be created
 */
bool mram_save(mram_info_t *mram, unsigned int dpu_id);

/**
 * @brief Loads an MRAM file into an MRAM image.
 *
 * @param mram the MRAM image
 * @param dpu_id index of the DPU
 */
void mram_load(mram_info_t *mram, unsigned int dpu_id);

/**
 * @brief Loads an MRAM file info into an MRAM image.
 *
 * @param mram the MRAM image
 * @param dpu_id index of the DPU
 */
void mram_load_info(mram_info_t *mram, unsigned int dpu_id);

#endif /* __INTEGRATION_MDPU_H__ */
