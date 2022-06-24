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

#include "common.h"
#include "index.h"

size_t mram_load(uint8_t **mram, unsigned int dpu_id);

void init_vmis(unsigned int nb_dpu, distribute_index_t *table);
void free_vmis(unsigned int nb_dpu);
void write_vmi(unsigned int num_dpu, unsigned int num_ref, coords_and_nbr_t *coords_and_nbr);

#endif /* __INTEGRATION_MDPU_H__ */
