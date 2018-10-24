#ifndef __VARTREE_H__
#define __VARTREE_H__

#define MAX_SIZE_ALLELE 10

typedef struct variant {
        char* chr;
        int offset;
        int pos;
        char ref[MAX_SIZE_ALLELE];
        char alt[MAX_SIZE_ALLELE];
        int depth;
} variant_t;

typedef struct variant_tree {
        int pos;
        variant_t* vars;
        struct variant_tree* right;
        struct variant_tree* left;
        int height;
} variant_tree_t;

void insert_variants(variant_tree_t **variant_list, variant_t* var);
// TODO: if not used at the end, remove this function
void visu_variants(variant_tree_t **variant_list);
void free_variant_tree(variant_tree_t *variant_list);

#endif /* __VARTREE_H__ */
