/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __COMPARE_H__
#define __COMPARE_H__

#include <stdint.h>
#include "upvc.h"

/**
 * @brief Structure to store the backtrack to be used after DPD algorithm.
 *
 * @var type  Difference between two neighbour sequence.
 * @var ix    X axis index of the path.
 * @var jx    Y axis index of the path.
 */
typedef struct {
        int type;
        int ix;
        int jx;
} backtrack_t;

/**
 * @brief Compute the alignment distance by dynamical programming on the diagonals of the matrix.
 *
 * @param s1          A neighbour sequence of int8_t.
 * @param s2          A neighbour sequence of int8_t.
 * @param backtrack   Output backtrack use to minimize the distance.
 * @param reads_info  Information on the size of the seed and the neighbour.
 *
 * @return the distance between s1 and s2.
 */
int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, reads_info_t *reads_info);

#endif /* __COMPARE_H__ */
