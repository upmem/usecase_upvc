#ifndef __CODE_H__
#define __CODE_H__

#include <stdint.h>
#include "upvc.h"

/*
 * Compute the code of a seed of size SIZE_SEED (see upvc.h).
 */
int code_seed(int8_t *sequence);
/*
 * Compute the code of a neighbour (of a seed).
 * The size is found in the "reads_info" argument.
 */
void code_neighbour(int8_t *sequence, int8_t *code, reads_info_t *reads_info);
/*
 * Compute the sequence of a neighbour from its code.
 * The size is found in the "reads_info" argument.
 * TODO: if not used at the end, remove this function
 */
void decode_neighbour(int8_t *code, int8_t *sequence, reads_info_t *reads_info);

#endif /* __CODE_H__ */
