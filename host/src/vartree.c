/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upvc.h"
#include "vartree.h"

/**
 * @brief Most of this code is inspired by geeksforgeeks.org
 */

static int max(int a, int b) { return (a > b) ? a : b; }

static variant_tree_t *alloc_variant(int pos, variant_t *vars)
{
    variant_tree_t *variant_tree = malloc(sizeof(variant_tree_t));
    variant_tree->pos = pos;
    variant_tree->vars = vars;
    variant_tree->right = NULL;
    variant_tree->left = NULL;
    variant_tree->height = 1;
    return variant_tree;
}

static int get_height(variant_tree_t *variant_tree)
{
    if (variant_tree == NULL) {
        return 0;
    } else {
        return variant_tree->height;
    }
}

static int get_balance(variant_tree_t *variant_tree)
{
    if (variant_tree == NULL) {
        return 0;
    } else {
        return get_height(variant_tree->left) - get_height(variant_tree->right);
    }
}

/**
 * @brief A utility function to right rotate subtree rooted with variant_tree
 * See the diagram given above.
 */
static variant_tree_t *right_rotate(variant_tree_t *variant_tree)
{
    variant_tree_t *right_tree = variant_tree->left;
    variant_tree_t *left_tree = right_tree->right;

    /* Perform rotation */
    right_tree->right = variant_tree;
    variant_tree->left = left_tree;

    /* Update heights */
    variant_tree->height = max(get_height(variant_tree->left), get_height(variant_tree->right)) + 1;
    right_tree->height = max(get_height(right_tree->left), get_height(right_tree->right)) + 1;

    /* Return new root */
    return right_tree;
}

/**
 * @brief A utility function to left rotate subtree rooted with variant_tree
 * See the diagram given above.
 */
static variant_tree_t *left_rotate(variant_tree_t *variant_tree)
{
    variant_tree_t *right_tree = variant_tree->right;
    variant_tree_t *left_tree = right_tree->left;

    /* Perform rotation */
    right_tree->left = variant_tree;
    variant_tree->right = left_tree;

    /* Update heights */
    variant_tree->height = max(get_height(variant_tree->left), get_height(variant_tree->right)) + 1;
    right_tree->height = max(get_height(right_tree->left), get_height(right_tree->right)) + 1;

    /* Return new root */
    return right_tree;
}

static variant_tree_t *insert(variant_tree_t *variant_tree, variant_t *var)
{
    /* if the tree does not exist, initialize it with the given variant */
    if (variant_tree == NULL) {
        return alloc_variant(var->pos, var);
    }

    /* if the tree is rooted at a variant with the same position, increase the variant depth */
    if (variant_tree->pos == var->pos) {
        variant_tree->vars->depth += 1;
        free(var);
        return variant_tree;
    }

    if (var->pos > variant_tree->pos) { /* if the variant is at a higher position than the root,
                                         *then insert it in the right subtree */
        variant_tree->right = insert(variant_tree->right, var);
    } else { /* if the variant is at a lower position than the variant_tree, then insert it in the left subtree */
        variant_tree->left = insert(variant_tree->left, var);
    }

    variant_tree->height = 1 + max(get_height(variant_tree->left), get_height(variant_tree->right));

    int balance = get_balance(variant_tree);

    /* If balance magnitude is greater than 1, then the tree needs to be rebalanced */

    /* Left Left Case */
    if (balance > 1 && var->pos < variant_tree->left->pos) {
        return right_rotate(variant_tree);
    }

    /* Right Right Case */
    if (balance < -1 && var->pos > variant_tree->right->pos) {
        return left_rotate(variant_tree);
    }

    /* Left Right Case */
    if (balance > 1 && var->pos > variant_tree->left->pos) {
        variant_tree->left = left_rotate(variant_tree->left);
        return right_rotate(variant_tree);
    }

    /* Right Left Case */
    if (balance < -1 && var->pos < variant_tree->right->pos) {
        variant_tree->right = right_rotate(variant_tree->right);
        return left_rotate(variant_tree);
    }

    /* Tree is still balanced */
    return variant_tree;
}

void insert_variants(variant_tree_t **variant_list, variant_t *var) { *variant_list = insert(*variant_list, var); }

void free_variant_tree(variant_tree_t *variant_tree)
{
    if (variant_tree == NULL) {
        return;
    }
    free_variant_tree(variant_tree->left);
    free_variant_tree(variant_tree->right);
    free(variant_tree->vars);
    free(variant_tree);
}
