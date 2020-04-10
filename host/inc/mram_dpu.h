/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_MDPU_H__
#define __INTEGRATION_MDPU_H__

/**
 * @brief Mapping of the DPU MRAMs.
 *
 * Defines the structures representing the DPU MRAMs on both the host and DPU side.
 */
#include <stddef.h>
#include <stdint.h>

#define FILE_NAME_SIZE 24

/**
 * @brief create the name of the 'dpu_id'th mram
 *
 * @param dpu_id index of the mram
 *
 * @return the name of the mram
 */
char *make_mram_file_name(char *str, unsigned int dpu_id);

/**
 * @brief Loads an MRAM file into an MRAM image.
 *
 * @param mram the MRAM image
 * @param dpu_id index of the DPU
 *
 * @return the size of the mram
 */
size_t mram_load(uint8_t *mram, unsigned int dpu_id);

#endif /* __INTEGRATION_MDPU_H__ */
