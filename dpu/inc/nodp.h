/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_NODP_H__
#define __INTEGRATION_NODP_H__

#include <stdint.h>

/**
 * @brief Compare two sets of symbols to produce a score.
 *
 * This version is a fast comparison of two sequences.
 */

#define COST_SUB (10)

static const int TT[256] = {0, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40,
                            20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40};

/**
 * @brief Compares two sequences of symbols to assign a score, based on the number of substitutions and INDELs
 *
 * The function computes the distance between two sequences, if there are only substitutions, or the detection
 * of INDELs.
 *
 * @param s1           The first sequence, packed in a byte stream
 * @param s2           The second sequence, packed in a byte stream
 * @param nbr_len      The size of each vector, in bytes
 * @param max_score    A threshold, beyond which there is no need to go further (because we already have a better score)
 * @return -1 if INDELs are detected, otherwise a score specifying the distance between s1 and s2 with substitutions
 */
int noDP(uint8_t *s1, uint8_t *s2, unsigned int nbr_len, int max_score) {
        int i, j, x, v, V, V1, V2;
        int score = 0;
        for (i = 0; i < nbr_len; i++) {
                x = (int) (s1[i] ^ s2[i]);
                x = x & 0xFF;
                v = TT[x];
                if (v > COST_SUB) {/* More than one difference found */
                        j = i + 1;
                        if (j < nbr_len - 3) {
                                V1 = s1[j] | (s1[j + 1] << 8) | (s1[j + 2] << 16) | (s1[j + 3] << 24);
                                V2 = s2[j] | (s2[j + 1] << 8) | (s2[j + 2] << 16) | (s2[j + 3] << 24);
                                /* V = (V1 ^ V2) & 0xFFFFFF; /\* Check if the next 8 characters are the same *\/ */
                                /* if (V != 0) {              /\* No, check INDELs                            *\/ */
                                V = (V1 ^ (V2 >> 2)) & 0x3FFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V1 ^ (V2 >> 4)) & 0xFFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V1 ^ (V2 >> 6)) & 0x3FFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V1 ^ (V2 >> 8)) & 0xFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }

                                V = (V2 ^ (V1 >> 2)) & 0x3FFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V2 ^ (V1 >> 4)) & 0xFFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V2 ^ (V1 >> 6)) & 0x3FFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                                V = (V2 ^ (V1 >> 8)) & 0xFFFFFF;
                                if (V == 0) {
                                        return -1;
                                }
                        }
                }

                score += v;
                if (score > max_score) {
                        break;
                }
        }
        return score;
}

#endif //INTEGRATION_NODP_H
