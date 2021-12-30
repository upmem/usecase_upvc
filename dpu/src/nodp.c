/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdint.h>

#include "common.h"
#include <defs.h>

#define COST_SUB 10

const uint8_t translation_table[256] = { 0, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20, 10, 20, 20, 20, 20, 30,
    30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20,
    30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
    30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30,
    30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30,
    30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20,
    30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 10, 20, 20, 20, 20, 30, 30, 30, 20, 30, 30, 30, 20, 30, 30, 30,
    20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40,
    40, 20, 30, 30, 30, 30, 40, 40, 40, 30, 40, 40, 40, 30, 40, 40, 40 };

#define CMP(a, b, shift) {                                                      \\
	if (a>>shift == (b & (0xffffffff>>shift))) {                            \\
		return UINT_MAX;                                                \\
	}

#define CMP_PAIR(a,b,shift) {                                                   \\
	CMP(a,b,shift)                                                          \\
	CMP(b,a,shift)
#if 0


uint32_t nodp(uint8_t *s1, uint8_t *s2, uint32_t max_score, uint32_t size_neighbour_in_bytes)
{
    uint32_t score = 0;
    uint32_t i = 0;
    while (i < size_neighbour_in_bytes) {
        uint32_t s1_word = s1[0];
        uint32_t s2_word = s2[0];
	s1+=sizeof(uint32_t);
	s2+=sizeof(uint32_t);
        uint32_t s_xor = s1_word ^ s2_word;
        for (uint32_t k = 0; (k < sizeof(uint32_t)) && (i < size_neighbour_in_bytes); k++, i++) {
            uint32_t s_translated = translation_table[s_xor & 0xFF];
            s_xor >>= CHAR_BIT;
        }
    }
    if (score>MAX_SCORE) {
	    uint32_t last_word1 = s1[0];
	    uint32_t last_word2 = s2[0];
	    CMP_PAIR(last_word1, last_word2, 2)
	    CMP_PAIR(last_word1, last_word2, 4)
	    CMP_PAIR(last_word1, last_word2, 6)
	    CMP_PAIR(last_word1, last_word2, 8)
    }
    return score;
}

#endif
