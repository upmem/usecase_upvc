/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __VARTREE_H__
#define __VARTREE_H__

#define MAX_SIZE_ALLELE 8

#include <stdint.h>

typedef struct variant {
    uint32_t score;
    uint32_t depth;
    char ref[MAX_SIZE_ALLELE];
    char alt[MAX_SIZE_ALLELE];
    int8_t *reads;
    struct variant *next;
} variant_t;

void variant_tree_insert(variant_t *var, uint32_t seq_nr, uint32_t offset_in_chr, int8_t *read);

void variant_tree_init();
void variant_tree_free();

void create_vcf();

#endif /* __VARTREE_H__ */
