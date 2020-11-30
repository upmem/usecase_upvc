/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <attributes.h>

#include "common.h"
#include "odpd.h"

/* Need big space to store one triplet of matrices per tasklet */
__dma_aligned uint8_t __M[3 * NR_TASKLETS * SIZEOF_MATRIX(NB_BYTES_TO_SYMS(SIZE_NEIGHBOUR_IN_BYTES, 0))];
