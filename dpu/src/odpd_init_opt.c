/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "odpd.h"
#include "common.h"

/* Need big space to store one triplet of matrices per tasklet */
uint8_t __M[3 * NB_TASKLET_PER_DPU * SIZEOF_MATRIX(NB_BYTES_TO_SYMS(SIZE_NEIGHBOUR_IN_BYTES, 0))];
