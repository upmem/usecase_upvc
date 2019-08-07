/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "odpd.h"
#include "debug.h"
#include <alloc.h>

#define NB_DIAGS 15
#define COST_SUB 10
#define COST_GAPO 11
#define COST_GAPE 1

static inline int min(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}

/* Need three huge tables per tasklet to work. */
static int *__D;
static int *__P;
static int *__Q;

static inline unsigned int nb_items_per_matrix(unsigned int nbr_sym_len)
{
    /* One matrix is ...[2][nbr_sym_len + 1] */
    return (nbr_sym_len + 1) << 1;
}

static inline unsigned int sizeof_matrix(unsigned int nbr_sym_len)
{
    /* Each matrix contains ints = 4 bytes */
    return nb_items_per_matrix(nbr_sym_len) << 2;
}

void odpd_init(unsigned int nb_tasklets, unsigned int nbr_sym_len)
{
    __D = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
    __P = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
    __Q = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
}

static inline int *get_matrix_for_tasklet(int *base, unsigned int tid, unsigned int nbr_sym_len)
{
    uint8_t *as_bytes = (uint8_t *)base;
    uint8_t *result = as_bytes + sizeof_matrix(nbr_sym_len) * tid;
    return (int *)result;
}

static inline int *_at(int *M, int x, int y, int nbr_sym_len)
{
    /* x is equal to 0 or 1 */
    return M + (((x == 0) ? 0 : nbr_sym_len + 1) + y);
}

#define D(x, y) *_at(_D, x, y, nbr_sym_len)
#define P(x, y) *_at(_P, x, y, nbr_sym_len)
#define Q(x, y) *_at(_Q, x, y, nbr_sym_len)

int odpd(const uint8_t *s1, const uint8_t *s2, int max_score, unsigned int nbr_sym_len)
{
    int *_D = get_matrix_for_tasklet(__D, tid, nbr_sym_len);
    int *_P = get_matrix_for_tasklet(__P, tid, nbr_sym_len);
    int *_Q = get_matrix_for_tasklet(__Q, tid, nbr_sym_len);
    unsigned int tid = me();

    int i, j, d, lp, pp, QP, min_score;

    for (j = 0; j <= NB_DIAGS / 2 + 1; j++) {
        P(0, j) = 99;
        Q(0, j) = 99;
    }
    P(1, 0) = 99;
    Q(1, 0) = 99;
    for (j = 0; j <= NB_DIAGS / 2 + 1; j++) {
        D(0, j) = j * COST_SUB;
    }

    for (i = 1; i < NB_DIAGS / 2 + 1; i++) {
        min_score = 99;
        pp = i & 1; /* i % 2 */
        lp = pp ^ 1; /* (i - 1) % 2 */
        D(pp, 0) = i * COST_SUB;
        for (j = 1; j <= i + NB_DIAGS / 2; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            if (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) != ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3)) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);
            if (D(pp, j) < min_score) {
                min_score = D(pp, j);
            }
        }
        Q(pp, j) = 99;
        D(pp, j) = 99;
        if (min_score > max_score) {
            return min_score;
        }
    }
    for (i = NB_DIAGS / 2 + 1; i < nbr_sym_len - NB_DIAGS / 2; i++) {
        min_score = 99;
        pp = i & 1; /* i % 2 */
        lp = pp ^ 1; /* (i - 1) % 2 */
        j = i - NB_DIAGS / 2 - 1;
        P(pp, j) = 99;
        D(pp, j) = 99;
        for (j = i - NB_DIAGS / 2; j <= i + NB_DIAGS / 2; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            if (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) != ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3)) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);
            if (D(pp, j) < min_score) {
                min_score = D(pp, j);
            }
        }
        Q(pp, j) = 99;
        D(pp, j) = 99;
        if (min_score > max_score) {
            return min_score;
        }
    }
    min_score = 99;

    for (i = nbr_sym_len - NB_DIAGS / 2; i < nbr_sym_len + 1; i++) {
        pp = i & 1; /* i % 2 */
        lp = pp ^ 1; /* (i - 1) % 2 */
        j = i - NB_DIAGS / 2 - 1;
        P(pp, j) = 99;
        D(pp, j) = 99;
        for (j = i - NB_DIAGS / 2; j < nbr_sym_len + 1; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            if (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) != ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3)) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);
        }
        if (D(pp, nbr_sym_len) < min_score) {
            min_score = D(pp, nbr_sym_len);
        }
    }
    i = nbr_sym_len;
    pp = i & 1; /* i % 2 */
    for (j = i - NB_DIAGS / 2; j < nbr_sym_len + 1; j++) {
        if (D(pp, j) < min_score) {
            min_score = D(pp, j);
        }
    }
    return min_score;
}
