/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "upvc.h"
#include "vartree.h"
#include "genome.h"

static variant_t **variant_list[MAX_SEQ_GEN] = {NULL};
static pthread_mutex_t mutex;

void variant_tree_init() {
    genome_t *genome = genome_get();
    pthread_mutex_init(&mutex, NULL);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        variant_list[each_seq] = (variant_t **)calloc(genome->len_seq[each_seq], sizeof(variant_t *));
    }
}

void variant_tree_insert(variant_t *var, uint32_t seq_nr, uint32_t offset_in_chr) {
    pthread_mutex_lock(&mutex);
    variant_t **entry = &variant_list[seq_nr][offset_in_chr];
    variant_t *vars = *entry;
    while (vars != NULL) {
        if (!strcmp(vars->ref, var->ref) && !strcmp(vars->alt, var->alt)) {
            vars->depth++;
            vars->score += var->score;
            free(var);
            goto end;
        }
        vars = vars->next;
    }
    var->next = *entry;
    *entry = var;

end:
    pthread_mutex_unlock(&mutex);
}

void variant_tree_free() {
    genome_t *genome = genome_get();
    pthread_mutex_destroy(&mutex);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        for (unsigned int i = 0; i < genome->len_seq[each_seq]; i++) {
            variant_t *tmp = variant_list[each_seq][i];
            while (tmp != NULL) {
                variant_t *to_free = tmp;
                tmp = tmp->next;
                free(to_free);
            }
        }
        free(variant_list[each_seq]);
    }
}

variant_t ***variant_tree_get() { return variant_list; }
