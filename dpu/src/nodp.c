#include <stdint.h>

/**
 * @brief Compare two sets of symbols to produce a score.
 *
 * This version is a fast comparison of two sequences.
 */

#define COST_SUB (10)

static const int translation_table[256] = {0, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20,
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

int noDP(uint8_t *s1, uint8_t *s2, unsigned int nbr_len, unsigned int delta, int max_score) {
        int score = 0;
        for (int i = 0; i < nbr_len - delta; i++) {
                int s_xor = (int) (s1[i] ^ s2[i]) & 0xFF;
                int s_translated = translation_table[s_xor];
                if (s_translated > COST_SUB) {/* More than one difference found */
                        int j = i + 1;
                        if (j < nbr_len - delta - 3) {
                                int V1 = s1[j] | (s1[j + 1] << 8) | (s1[j + 2] << 16) | (s1[j + 3] << 24);
                                int V2 = s2[j] | (s2[j + 1] << 8) | (s2[j + 2] << 16) | (s2[j + 3] << 24);
                                int V = (V1 ^ (V2 >> 2)) & 0x3FFFFFFF;
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

                score += s_translated;
                if (score > max_score) {
                        break;
                }
        }
        return score;
}
