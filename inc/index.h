#ifndef __INDEX_H__
#define __INDEX_H__

#include "genome.h"
#include "upvc.h"

#define NB_SEED (16777216)

typedef struct index_seed {
        int nb_nbr;              /* Number of neighbour                                              */
        int offset;              /* Address in the DPU memory of the first neighbour to compute      */
        int num_dpu;             /* DPU number where the reference seed that match has been dispatch */
        struct index_seed *next; /* Needed to have a linked-list of seed informations                */
} index_seed_t;

/*
 * Create the linked-list of seed, and write them in the DPUs memories.⎈
 */
index_seed_t **index_genome(genome_t *ref_genome, times_ctx_t *times_ctx, reads_info_t *reads_info);
/*
 * Free the linked-list of seed.
 */
void free_index(index_seed_t **index_seed);

#endif /* __INDEX_H__ */
