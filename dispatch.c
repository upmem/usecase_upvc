#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

void dispatch_read(INDEX_SEED **SEED, int8_t* BUF_READ, int nb_read, TIMES *CT)
{
  int i, k, v, numdpu, num_read;
  int8_t *read;
  double t1, t2;
  INDEX_SEED *C;
  int8_t *buf = (int8_t *) malloc(sizeof(int8_t)*SIZE_NBR);
  int COUNT_READ[NB_DPU];

  t1 = my_clock();
  for (i=0; i<NB_DPU; i++) COUNT_READ[i] = 0;
  for (num_read=0; num_read<nb_read; num_read++)
    {
      read = &BUF_READ[num_read*SIZE_READ];
      // dispatch du read dans les DPUs
      v = code_seed(read);
      C = SEED[v];
      while (C!=NULL)
	{
	  k = COUNT_READ[C->num_dpu]; // k = nombre de reads deja memorise dans le DPU
	  // a remplacer par les transactions memoire correspondantes
	  write_count(C->num_dpu,k,C->nb_nbr);
	  write_offset(C->num_dpu,k,C->offset);
	  write_num(C->num_dpu,k,num_read); 
	  code_neighbor(&read[SIZE_SEED],buf);
	  write_neighbor_read(C->num_dpu,k,buf);
	  COUNT_READ[C->num_dpu]++;
	  if (COUNT_READ[C->num_dpu] >= MAX_NB_DPU_READ-1) { printf ("\nBuffer full (DPU %d)\n",C->num_dpu); exit(255); }
	  C = C->next;
	}
    }
  for (numdpu=0; numdpu < NB_DPU; numdpu++)
    {
      k = COUNT_READ[numdpu];
      write_num(numdpu,k,-1);
    }

  free(buf);
  t2 = my_clock();
  CT->dispatch_read = t2-t1;
  CT->tot_dispatch_read += t2-t1;
}
