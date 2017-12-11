#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"


GENOME *get_genome(char* name, TIMES *CT)
{
  
  FILE *fgen;          // fichier genome
  long sizefile;       // taille du fichier
  int  i, j, k;
  long x;
  double t1, t2;
  char BUF[1024];
  char filename[1024];
  GENOME *G;

  // lecture du fichier qui contient le genome (format FASTA)
  // ecriture des sequences dans GENOME

  t1 = my_clock();
  sprintf(filename,"%s.fasta",name);
  fgen = fopen(filename,"r");
  fseek(fgen,0,SEEK_END);
  sizefile = ftell(fgen);
  rewind(fgen);

  G = (GENOME *) malloc(sizeof(GENOME));
  G->sizefile = sizefile;
  G->data    = (int8_t *) malloc(sizeof(int8_t) * sizefile); 
  G->pt_seq  = (int *) malloc(sizeof(int) * MAX_SEQ_GEN);
  G->len_seq = (int *) malloc(sizeof(int) * MAX_SEQ_GEN);
  G->nb_seq = 0;
  x = 0;
  while (fgets(BUF,1024,fgen)!=NULL)
    {
      if (BUF[0] == '>')
	{
	  G->pt_seq[G->nb_seq] = x;
	  G->len_seq[G->nb_seq] = 0;
	  G->nb_seq++;
	}
      else
	{
	  k = strlen(BUF);
	  for (i=0; i<k-1; i++) 
	    {
	      j = (int) BUF[i];
	      G->data[x++] = (j>>1)&3; // A -> 0, C -> 1, G -> 3, T -> 2
	      G->len_seq[G->nb_seq-1]++;
	      }
	  }
      }
  fclose(fgen);
  t2 = my_clock();
  CT->get_genome = t2-t1;
  return G;
}

void free_genome(GENOME *G)
{
  free(G->pt_seq); 
  free(G->len_seq); 
  free(G->data); 
  free(G);
}
