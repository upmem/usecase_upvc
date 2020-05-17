/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_NODP_H__
#define __INTEGRATION_NODP_H__

#include <stdint.h>

/**
 * @brief Compares two sequences of symbols to assign a score, based on the number of substitutions and INDELs
 *
 * The function computes the distance between two sequences, if there are only substitutions, or the detection
 * of INDELs.
 *
 * @return -1 if INDELs are detected, otherwise a score specifying the distance between s1 and s2 with substitutions
 */
uint32_t nodp(uint8_t *s1, uint8_t *s2, uint32_t max_score, uint32_t size_neighbour_in_bytes);

#endif /* __INTEGRATION_NODP_H__ */
