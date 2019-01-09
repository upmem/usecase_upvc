/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdint.h>
#include "upvc.h"

#define CODE_SIZE (4)
#define MASK (CODE_SIZE - 1)

int code_seed(int8_t *sequence)
{
        int seed = 0;
        for (int i = 0; i < SIZE_SEED; i++) {
                if (sequence[i] >= CODE_SIZE) {
                        return -1;
                }
                seed = (seed * CODE_SIZE) + sequence[i];
        }
        return seed;
}

void code_neighbour(int8_t *sequence, int8_t *code, reads_info_t *reads_info)
{
        for (int i = 0; i < reads_info->size_neighbour_in_bytes; i++) {
                int j = i * 4;
                code[i] = ((sequence[j + 3] & MASK) * (CODE_SIZE * CODE_SIZE * CODE_SIZE))
                        + ((sequence[j + 2] & MASK) * (CODE_SIZE * CODE_SIZE))
                        + ((sequence[j + 1] & MASK) * (CODE_SIZE))
                        + ((sequence[j    ] & MASK));
        }
}
