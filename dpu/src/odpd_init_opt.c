#include <alloc.h>
#include "odpd.h"

/* Need big space to store one triplet of matrices per tasklet */
int *__M = 0;

void odpd_init(unsigned int nb_tasklets, unsigned int nbr_sym_len)
{
        unsigned int toto = 3 * SIZEOF_MATRIX(nbr_sym_len) * nb_tasklets;
        __M = mem_alloc(3 * SIZEOF_MATRIX(nbr_sym_len) * nb_tasklets);
}
