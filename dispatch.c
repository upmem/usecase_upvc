#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

void writeMemDPU(INDEX_SEED *C, int *COUNT_READ, int8_t *read, int num_read)
{
  int k;
  int8_t *buf = (int8_t *) malloc(sizeof(int8_t)*SIZE_NBR);

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
  free(buf);
}

void dispatch_read(INDEX_SEED **SEED, int8_t* BUF_READ, int nb_read, TIMES *CT)
{
  int i, k, v0, v1, v2, v3, numdpu, num_read;
  int8_t *read0, *read1, *read2, *read3;
  double t1, t2;
  INDEX_SEED *C0, *C1, *C2, *C3;
  int COUNT_READ[NB_DPU];

  t1 = my_clock();
  for (i=0; i<NB_DPU; i++) COUNT_READ[i] = 0;
  for (num_read=0; num_read<nb_read; num_read=num_read+4)
    {
      read0 = &BUF_READ[num_read*SIZE_READ];
      v0 = code_seed(read0);
      C0 = SEED[v0];

      read1 = &BUF_READ[(num_read+1)*SIZE_READ];
      v1 = code_seed(read1);
      C1 = SEED[v1];

      read2 = &BUF_READ[(num_read+2)*SIZE_READ];
      v2 = code_seed(read2);
      C2 = SEED[v2];

      read3 = &BUF_READ[(num_read+3)*SIZE_READ];
      v3 = code_seed(read3);
      C3 = SEED[v3];

      writeMemDPU(C0,COUNT_READ,read0,num_read);
      writeMemDPU(C2,COUNT_READ,read2,num_read+2);

      writeMemDPU(C1,COUNT_READ,read1,num_read+1);
      writeMemDPU(C3,COUNT_READ,read3,num_read+3);
    }
  for (numdpu=0; numdpu < NB_DPU; numdpu++)
    {
      k = COUNT_READ[numdpu];
      write_num(numdpu,k,-1);
    }

  t2 = my_clock();
  CT->dispatch_read = t2-t1;
  CT->tot_dispatch_read += t2-t1;
}
