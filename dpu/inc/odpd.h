/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_ODPD_H__
#define __INTEGRATION_ODPD_H__

#include <stdint.h>

/**
 * @brief Compare two sets of symbols to produce a score.
 *
 * The module implements a version of Smith&Waterman algorithm. It first requires initializing
 * memory resources, before invoking the comparison function odpd.
 */

/**
 * @brief Sets up an operating environment for the comparator to run on several tasklets in parallel.
 *
 * @param nb_tasklets  How many tasklets will work.
 * @param nbr_sym_len  Size of a neighbour, in number of symbols.
 */
void odpd_init(unsigned int nb_tasklets, unsigned int nbr_sym_len);

/**
 * @brief Compares two sequences of symbols to assign a score.
 *
 * Optimal Dynamic Programming Diagonal. The sequences are expressed a byte streams, i.e. 4
 * nucleotides per bytes.
 *
 * @param s1           The first vector
 * @param s2           The second vector
 * @param max_score    Any score above this threshold is good
 * @param nbr_sym_len  The number of symbols in s1 and s2
 * @param tid          Sysname of the invoking tasklet
 * @return A score
 */
int odpd(const uint8_t *s1, const uint8_t *s2, int max_score, unsigned int nbr_sym_len, unsigned int tid);

/* Assembler version */
int odpdx(const uint8_t *s1, const uint8_t *s2, int max_score, unsigned int nbr_sym_len, unsigned int tid);
#endif /* __INTEGRATION_ODPD_H__ */
