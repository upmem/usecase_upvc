/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INDEX_H__
#define __INDEX_H__

#include <stdint.h>

/**
 * @brief Structure of a list of index of neighbour that shared the same seed.
 *
 * @var nb_nbr   Number of neighbour.
 * @var offset   Address in the DPU memory of the first neighbour to compute.
 * @var num_dpu  DPU number where the reference seed that match has been dispatch.
 * @var next     Needed to have a linked-list of index.
 */
typedef struct index_seed {
    uint32_t nb_nbr;
    uint32_t offset;
    uint32_t num_dpu;
    struct index_seed *next;
} index_seed_t;

void index_save();

void index_load();

void index_init();

void index_free();

index_seed_t *index_get(int8_t *read);

unsigned int index_get_nb_dpu();

void index_copy_neighbour(int8_t *dst, int8_t *src);

#endif /* __INDEX_H__ */
