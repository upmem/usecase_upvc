/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <alloc.h>
#include "odpd.h"

/* Need big space to store one triplet of matrices per tasklet */
int *__M = 0;

void odpd_init(unsigned int nb_tasklets, unsigned int nbr_sym_len)
{
        __M = mem_alloc(3 * SIZEOF_MATRIX(nbr_sym_len) * nb_tasklets);
}
