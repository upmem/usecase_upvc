/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __VARTREE_H__
#define __VARTREE_H__

#define MAX_SIZE_ALLELE 10

#include <stdint.h>

typedef struct variant {
    char *chr;
    uint64_t offset;
    int64_t pos;
    char ref[MAX_SIZE_ALLELE];
    char alt[MAX_SIZE_ALLELE];
    int depth;
} variant_t;

typedef struct variant_tree {
    int64_t pos;
    variant_t *vars;
    struct variant_tree *right;
    struct variant_tree *left;
    int height;
} variant_tree_t;

void insert_variants(variant_t *var);

void free_variant_tree();

variant_tree_t *variant_list_get();

#endif /* __VARTREE_H__ */
