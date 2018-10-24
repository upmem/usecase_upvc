#ifndef __COMPARE_H__
#define __COMPARE_H__

#include <stdint.h>
#include "upvc.h"

typedef struct {
        int type;
        int ix;
        int jx;
} backtrack_t;

// TODO: if not used at the end, remove this function
void display_DPD_matrix(int **M, reads_info_t *reads_info);
/*
 * Compute the alignment distance by dynamical programming on the diagonals of the matrix.
 * Return the backtrack.
 */
int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, reads_info_t *reads_info);

#endif /* __COMPARE_H__ */
