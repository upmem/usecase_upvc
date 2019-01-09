/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <alloc.h>
#include <ktrace.h>
#include "odpd.h"

#define NB_DIAGS         15
#define COST_SUB         10
#define COST_GAPO        11
#define COST_GAPE        1
#define COST_INIT        99

extern int *__M;

static inline int min(int a, int b) { if (a < b) return a; else return b; }

static inline int *get_matrix_for_tasklet(unsigned int tid, unsigned int nbr_sym_len)
{
        return (int *) ((uint8_t *) __M + ( 3 * SIZEOF_MATRIX(nbr_sym_len) ) * tid);
}

#define LINE_SIZE 6
#define d0off  0
#define p0off  2
#define p0off1 3
#define q0off  4
#define q0off1 5
#define d1off  6
#define p1off  8
#define q1off  10

int odpd(const uint8_t *s1, const uint8_t *s2, int max_score, unsigned int nbr_sym_len, unsigned int tid)
{
        /* Loop variables */
        int *_M = get_matrix_for_tasklet(tid, nbr_sym_len);
        int *_Mpp, *_Mlp;
        unsigned int i, j, min_score;
        /* local variables */
        int cost;
        unsigned int d, QP;

        cost = 0;
        for (_Mpp = _M, _Mlp = _M + (NB_DIAGS / 2 + 1) * LINE_SIZE; _Mpp<=_Mlp; _Mpp+=LINE_SIZE) {
                _Mpp[d0off] = cost;
                _Mpp[p0off] = COST_INIT;
                _Mpp[q0off] = COST_INIT;
                cost += COST_SUB;
        }
        _M[p0off1] = COST_INIT;
        _M[q0off1] = COST_INIT;

        _Mpp = _M + 1;
        _Mlp = _M;

        /* Phase 1 */
        for (i = 1; i < NB_DIAGS / 2 + 1; i++) {
                uint8_t v1 = s1[(i - 1) / 4U] >> (2 * ((i - 1) % 4U));
                min_score = COST_INIT;
                _Mpp[d0off] = i*COST_SUB;
                for (j = 1; j <= i + NB_DIAGS / 2; j++,_Mpp+=LINE_SIZE,_Mlp+=LINE_SIZE) {
                        _Mpp[p1off] = min(_Mpp[d0off] + COST_GAPO, _Mpp[p0off] + COST_GAPE);
                        _Mpp[q1off] = min(_Mlp[d1off] + COST_GAPO, _Mlp[q1off] + COST_GAPE);
                        QP = min(_Mpp[p1off], _Mpp[q1off]);
                        d = _Mlp[d0off];
                        uint8_t v2 = s2[(j - 1) / 4U] >> (2 * ((j - 1) % 4U));
                        if ( ((v1^v2) & 3) != 0 ) {
                                d += COST_SUB;
                        }
                        _Mpp[d1off] = min(d, QP);
                        if (_Mpp[d1off] < min_score) {
                                min_score = _Mpp[d1off];
                        }
                }
                _Mpp[q1off] = COST_INIT;
                _Mpp[d1off] = COST_INIT;
                if (min_score > max_score) {
                        return min_score;
                }
                _Mpp = _M + ((i&1)^1);
                _Mlp = _M + ((i&1));
        }

        /* Phase 2 */
        for (; i < nbr_sym_len - NB_DIAGS / 2; i++) {
                uint8_t v1 = s1[(i - 1) / 4U] >> (2 * ((i - 1) % 4U));
                min_score = COST_INIT;
                _Mpp[p0off] = COST_INIT;
                _Mpp[d0off] = COST_INIT;
                for (j = i - NB_DIAGS / 2; j <= i + NB_DIAGS / 2; j++,_Mpp+=LINE_SIZE,_Mlp+=LINE_SIZE) {
                        _Mpp[p1off] = min(_Mpp[d0off] + COST_GAPO, _Mpp[p0off] + COST_GAPE);
                        _Mpp[q1off] = min(_Mlp[d1off] + COST_GAPO, _Mlp[q1off] + COST_GAPE);
                        QP = min(_Mpp[p1off], _Mpp[q1off]);
                        d = _Mlp[d0off];
                        uint8_t v2 = s2[(j - 1) / 4U] >> (2 * ((j - 1) % 4U));
                        if ( ((v1^v2) & 3) != 0 ) {
                                d += COST_SUB;
                        }
                        _Mpp[d1off] = min(d, QP);
                        if (_Mpp[d1off] < min_score) {
                                min_score = _Mpp[d1off];
                        }
                }
                _Mpp[q1off] = COST_INIT;
                _Mpp[d1off] = COST_INIT;
                _Mpp -= LINE_SIZE * 2 * (NB_DIAGS/2);
                _Mlp -= LINE_SIZE * 2 * (NB_DIAGS/2);
                _Mpp = (int*)((unsigned int)_Mpp^4);
                _Mlp = (int*)((unsigned int)_Mlp^4);
                if (min_score > max_score) {
                        return min_score;
                }
        }

        /* Phase 3 */
        min_score = COST_INIT;
        _M = _Mpp+LINE_SIZE;
        for (; i < nbr_sym_len + 1; i++) {
                uint8_t v1 = s1[(i - 1) / 4U] >> (2 * ((i - 1) % 4U));
                _Mpp[p0off] = COST_INIT;
                _Mpp[d0off] = COST_INIT;
                for (j = i - NB_DIAGS / 2; j < nbr_sym_len + 1; j++,_Mpp+=LINE_SIZE,_Mlp+=LINE_SIZE) {
                        _Mpp[p1off] = min(_Mpp[d0off] + COST_GAPO, _Mpp[p0off] + COST_GAPE);
                        _Mpp[q1off] = min(_Mlp[d1off] + COST_GAPO, _Mlp[q1off] + COST_GAPE);
                        QP = min(_Mpp[p1off], _Mpp[q1off]);
                        d = _Mlp[d0off];
                        uint8_t v2 = s2[(j - 1) / 4U] >> (2 * ((j - 1) % 4U));
                        if ( ((v1^v2) & 3) != 0 ) {
                                d += COST_SUB;
                        }
                        _Mpp[d1off] = min(d, QP);
                }
                if (_Mpp[d0off] < min_score) {
                        min_score = _Mpp[d0off];
                }
                _Mlp -= (_Mpp-_M);
                _Mpp = _M;
                _Mpp = (int*)((unsigned int)_Mpp^4);
                _Mlp = (int*)((unsigned int)_Mlp^4);
                _M =_Mpp+LINE_SIZE;
        }

        _Mpp = (int*)((unsigned int)_Mpp^4);
        for (j = nbr_sym_len - NB_DIAGS / 2; j < nbr_sym_len + 1; j++, _Mpp+=LINE_SIZE) {
                if (_Mpp[d0off] < min_score) min_score = _Mpp[d0off];
        }

        return min_score;
}
