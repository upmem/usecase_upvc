/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INDEX_H__
#define __INDEX_H__

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/queue.h>

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

TAILQ_HEAD(distribute_index_list, distribute_index);
typedef struct distribute_index {
    uint64_t workload;
    uint32_t size;
    uint32_t dpu_id;
    TAILQ_ENTRY(distribute_index) entries;
} distribute_index_t;

char *get_index_folder();

void index_load();

void index_create();

void index_create_folder();

void index_free();

index_seed_t *index_get(int8_t *read);

unsigned int index_get_nb_dpu();

void index_copy_neighbour(int8_t *dst, int8_t *src);

enum xfer_direction {
    xfer_read,
    xfer_write,
};
static inline void xfer_file(uint8_t *buffer, size_t size_to_xfer, FILE *f, enum xfer_direction direction)
{
    uint64_t index = 0ULL;
    while (size_to_xfer != 0) {
        size_t size_xfer;
        if (direction == xfer_read) {
            size_xfer = fread(&buffer[index], 1, size_to_xfer, f);
        } else {
            size_xfer = fwrite(&buffer[index], 1, size_to_xfer, f);
        }
        assert(size_xfer > 0);
        size_to_xfer -= size_xfer;
        index += size_xfer;
    }
}

#endif /* __INDEX_H__ */
