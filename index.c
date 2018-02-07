#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

static int cmplong (void const *a, void const *b)
{
  long pa =  *(const long *) a;
  long pb =  *(const long *) b;
  if (pa<pb) return 1; else return -1;
}

INDEX_SEED **index_genome(GENOME *RG, TIMES *CT)
{
  double t1,t2;
  int seed, d,i, x, z, v, k, numdpu;
  long c, ns, lx, is;
  INDEX_SEED *C;
  long workload[NB_DPU];
  int sizeidx[NB_DPU];
  int OFFSET[NB_DPU];
  long *TMP;
  int *COUNT;
  int *CC;
  int8_t *buf;

  INDEX_SEED **SEED = (INDEX_SEED **) malloc(sizeof(INDEX_SEED *)*NB_SEED);

  t1 = my_clock();
 
  // comptage des graines
  COUNT  = (int *) malloc(sizeof(int)*NB_SEED);
  for (i=0; i<NB_SEED; i++) COUNT[i]=0; // initialisation de tous les compteurs a zero
  for (ns=0; ns<RG->nb_seq; ns++) // pour chaque sequence du genome
    {
      lx = RG->pt_seq[ns]; // lx = position dans GENOME du 1er caractere de la sequence numero ns
      for (is=0; is<RG->len_seq[ns]-SIZE_NBR-SIZE_SEED+1; is++) // pour chaque graine de la sequence
	{
	  v=code_seed(&RG->data[lx+is]); // v contient le code de la graine
	  if (v>=0) COUNT[v]++; // incrementation du compteur des graines
	}
    }

  // initialisation des index des graines a repartir sur les DPUs
  for (v=0; v<NB_SEED; v++)
    {
      x = (COUNT[v] / MAX_SIZE_IDX_SEED) + 1;             // pour chaque graine on calcule le nombre d'index
      k = (COUNT[v] / x)+1;                               // k = nombre de voisinage
      C = (INDEX_SEED *) malloc(sizeof(INDEX_SEED));      // on alloue en memoire le nombre d'index necessaire
      C->nb_nbr = k;
      C->next = NULL;
      SEED[v] = C;
      for (i=1; i<x; i++)
	{
	  C = (INDEX_SEED *) malloc(sizeof(INDEX_SEED));
	  C->nb_nbr = k;
	  C->next = SEED[v];
	  SEED[v] = C;
	}
      C->nb_nbr = COUNT[v]-((x-1)*k);
    }

  // repartition des index dans les DPUs

  // 1 - tri des graines par ordre decroissant
  TMP = (long *)malloc(sizeof(long)*NB_SEED);    // structure temporaire
  for (i=0; i<NB_SEED; i++) TMP[i] = (((long) COUNT[i]) << 32) + ((long) i); // aggregation count + seed
  qsort(TMP,NB_SEED,sizeof(long),cmplong);


  // 2 - repartition dans les DPUS
  for (numdpu=0; numdpu<NB_DPU; numdpu++) { workload[numdpu] = 0; sizeidx[numdpu] = 0; }

  d=0;
  for (i=0; i<NB_SEED; i++) 
    {
      seed = TMP[i] & 0xFFFFFFFF;  // on recupere les graines qui ont ete triees par ordre decroissant
      //printf ("%x %d\n",seed,COUNT[seed]);
      C = SEED[seed];
      while (C!= NULL)
	{
	  if (COUNT[seed]!=0)
	    {
	      c = 100000000000; k=0;
	      for (numdpu=0; numdpu<NB_DPU; numdpu++) 
		{
		  if (workload[d] < c) { k=d; c=workload[d]; }
		  d = (d+1)%NB_DPU;
		}
	    }
	  else
	    {
	      k=d;
	    }
	  sizeidx[k] += C->nb_nbr;
	  workload[k] += (long) (C->nb_nbr * COUNT[seed]);
	  C->num_dpu = k;
	  d = (d+1)%NB_DPU;
	  C = C->next;
	}
    }

  //data for ploting the size of the index
  //for (i=0; i<NB_DPU; i++) printf ("%d %ld %d\n",i,workload[i],sizeidx[i]*35);
  //exit (0);


  // calcul des offsets
  // initialisation du tableau IXH->OFFSET
  for (i=0; i<NB_DPU; i++) OFFSET[i]=0;        //pour chaque DPU on initialise l'offset a 0
  for (i=0; i<NB_SEED; i++)
    {
      C = SEED[i];
      while (C!=NULL)
	{
	  C->offset = OFFSET[C->num_dpu];      // mise a jour de l'offset
	  OFFSET[C->num_dpu] += C->nb_nbr;             // calcul de l'offset suivant
	  C = C->next;
	}
    }

  malloc_dpu();
  for (i=0; i<NB_DPU; i++) malloc_neighbor_idx(i,sizeidx[i]);


  // memorisation des sequences dans l'index des DPUs
  // le processeur écrit directement dans la memoire des DPUs
  CC = (int *) malloc(sizeof(int)*NB_SEED);
  buf = (int8_t *) malloc(sizeof(int8_t)*SIZE_NBR);
  for (i=0; i<NB_SEED; i++) CC[i] = 0;
  for (ns=0; ns<RG->nb_seq; ns++) // pour chaque sequence du genome
    {
      lx = RG->pt_seq[ns]; // x = position dans GENOME du 1er caractere de la sequence numero ns
      for (is=0; is<RG->len_seq[ns]-SIZE_NBR-SIZE_SEED+1; is++) // pour chaque graine de la sequence
	{
	  v=code_seed(&RG->data[lx+is]);
	  if (v<0) continue;
	  C = SEED[v];
	  z = 0;
	  while (C!= NULL)
	    {
	      if (CC[v] < C->nb_nbr + z) break;
	      z=z+C->nb_nbr;
	      C = C->next;
	    }
	  k = C->offset+CC[v]-z;
	  // codage et ecriture des voisinage dans le DPU concerne
	  // a remplacer par les transactions memoire correspondantes
	  code_neighbor(&RG->data[lx+is+SIZE_SEED],buf);
	  write_neighbor_idx(C->num_dpu,k,buf);
	  write_coordinate(C->num_dpu,k,((ns<<32)+is));
	  CC[v]++;
	}
    }

  free(TMP);
  free(COUNT);
  free(CC); 
  free(buf);
  t2 = my_clock();
  CT->index_genome = t2-t1;
  return SEED;
}

void free_index(INDEX_SEED **SEED)
{
  INDEX_SEED *C, *CC;
  int i;

  for (i=0; i<NB_SEED; i++)
    {
      C = SEED[i];
      while (C!=NULL)
	{
	  CC = C;
	  C = C->next;
	  free(CC);
	}
    }
  free(SEED);
}
