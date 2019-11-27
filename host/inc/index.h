/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INDEX_H__
#define __INDEX_H__

#define NB_SEED (1 << (SIZE_SEED << 1)) /* NB_SEED = 4 ^ (SIZE_SEED) */

#define SEED_FILE_LOG ("seeds.log")

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

/**
 * @brief Saves the index seeds into a text file, which can be reloaded for the stages follownig indexing.
 *
 * @param index_seed  Table of list of index to be save.
 */
void save_index_seeds(index_seed_t **index_seed);

/**
 * @brief Create a index seeds from a text file that has been created with "load_index_seeds".
 *
 * @return The table of list of index loaded from the text file.
 */
index_seed_t **load_index_seeds();

#include "dispatch.h"
#include "dpus_mgmt.h"
#include "genome.h"
#include "upvc.h"
#include "vmi.h"

#include "backends_functions.h"

/**
 * @brief Create the linked-list of index, and write them into the DPUs memories.
 *
 * @param ref_genome  Reference genome to use to create the index.
 * @param nb_dpu      Number of dpu available to compute.
 * @param times_ctx   Times information for the whole application.
 *
 * @return The Table of list of index created.
 */
index_seed_t **index_genome(genome_t *ref_genome, int nb_dpu, times_ctx_t *times_ctx);

/**
 * @brief Free the table of linked-list of index.
 *
 * @param index_seed The table to be freed.
 */
void free_index(index_seed_t **index_seed);

/**
 * @brief Print a table of list of index.
 *
 * @param index_seed  Table to be printed.
 * @param out         File in which to printed the table.
 */
void print_index_seeds(index_seed_t **SEED, FILE *out);

#endif /* __INDEX_H__ */
