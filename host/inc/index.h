#ifndef __INDEX_H__
#define __INDEX_H__

#include <stdbool.h>
#include "genome.h"
#include "upvc.h"

#define NB_SEED (16777216)

#define SEED_FILE_LOG ("seeds.log")

typedef struct index_seed {
        int nb_nbr;              /* Number of neighbour                                              */
        int offset;              /* Address in the DPU memory of the first neighbour to compute      */
        int num_dpu;             /* DPU number where the reference seed that match has been dispatch */
        struct index_seed *next; /* Needed to have a linked-list of seed information                 */
} index_seed_t;

/*
 * Saves the index seeds into a text file, which can be reloaded for the stages follownig indexing.
 */
void save_index_seeds(index_seed_t **index_seed);
/*
 * Create a index seeds from a text file that has been created with "load_index_seeds"
 */
index_seed_t **load_index_seeds();
/*
 * Create the linked-list of seed, and write them into the DPUs memories.
 */
index_seed_t **index_genome(genome_t *ref_genome,
                            int nb_dpu,
                            times_ctx_t *times_ctx,
                            reads_info_t *reads_info,
                            bool simulation_mode);
/*
 * Free the linked-list of seed.
 */
void free_index(index_seed_t **index_seed);
/*
 * print a index seed
 */
void print_index_seeds(index_seed_t **SEED, FILE *out, reads_info_t *reads_info);

#endif /* __INDEX_H__ */
