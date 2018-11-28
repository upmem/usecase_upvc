/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __CODE_H__
#define __CODE_H__

#include <stdint.h>
#include "upvc.h"

/**
 * @brief Compute the code of a seed of size SIZE_SEED.
 *
 * @param sequence  Sequence of int8_t of a seed.
 *
 * @return Code of the seed in a int.
 */
int code_seed(int8_t *sequence);

/**
 * @brief Compute the code of a neighbour (of a seed).
 *
 * @param sequence    Input neighbour as a table of int8_t.
 * @param code        Ouput sequence of int8_t of the neighbour.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
void code_neighbour(int8_t *sequence, int8_t *code, reads_info_t *reads_info);

#endif /* __CODE_H__ */