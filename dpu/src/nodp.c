/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdint.h>

#include <defs.h>
#include "common.h"

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

#if 0
#define CMP(V1, V2, shift)                                                                                                       \
    do {                                                                                                                         \
        if ((V2 >> shift) == (V1 & (0xffffffff >> shift))) {                                                                     \
            return UINT_MAX;                                                                                                     \
        }                                                                                                                        \
    } while (0)

#define CMP_PAIR(V1, V2, shift)                                                                                                  \
    do {                                                                                                                         \
        CMP(V1, V2, shift);                                                                                                      \
        CMP(V2, V1, shift);                                                                                                      \
    } while (0)

uint32_t nodp(uint8_t *s1, uint8_t *s2, uint32_t max_score, uint32_t size_neighbour_in_bytes)
{
    uint32_t score = 0;
    uint32_t s1_acc_h = *(uint32_t *)s1;
    uint32_t s2_acc_h = *(uint32_t *)s2;
    uint32_t i = 0;
    while (i < size_neighbour_in_bytes) {
        uint32_t s1_acc_l = s1_acc_h;
        uint32_t s2_acc_l = s2_acc_h;
        s1_acc_h = *(uint32_t *)(&s1[i + sizeof(uint32_t)]);
        s2_acc_h = *(uint32_t *)(&s2[i + sizeof(uint32_t)]);
        uint32_t s_xor = s1_acc_l ^ s2_acc_l;
        for (uint32_t k = 0; (k < sizeof(uint32_t)) && (i < size_neighbour_in_bytes); k++, i++) {
            uint32_t s_translated = translation_table[s_xor & 0xFF];
            if ((s_translated > COST_SUB) && ((i + sizeof(uint32_t)) < size_neighbour_in_bytes)) {
                uint32_t V1, V2;
                if (k != (sizeof(uint32_t) - 1)) {
                    uint32_t jr = (k + 1) * CHAR_BIT;
                    uint32_t jl = sizeof(uint32_t) * CHAR_BIT - jr;
                    V1 = (s1_acc_l >> jr) | ((s1_acc_h & ((1 << jr) - 1)) << jl);
                    V2 = (s2_acc_l >> jr) | ((s2_acc_h & ((1 << jr) - 1)) << jl);
                } else {
                    V1 = s1_acc_h;
                    V2 = s2_acc_h;
                }

                CMP_PAIR(V1, V2, 2);
                CMP_PAIR(V1, V2, 4);
                CMP_PAIR(V1, V2, 6);
                CMP_PAIR(V1, V2, 8);
            }

            score += s_translated;
            if (score > max_score) {
                return score;
            }
            s_xor >>= CHAR_BIT;
        }
    }
    return score;
}
#endif
