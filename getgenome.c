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
  G->pt_seq  = (long *) malloc(sizeof(long) * MAX_SEQ_GEN);
  G->len_seq = (int *) malloc(sizeof(int) * MAX_SEQ_GEN);
  G->name = (char **) malloc(sizeof(char*) * MAX_SEQ_GEN);
  G->nb_seq = 0;
  G->filename = strdup(filename);

  x = 0;
  while (fgets(BUF,1024,fgen)!=NULL)
    {
      if (BUF[0] == '>')
	{
	  G->pt_seq[G->nb_seq] = x;
	  G->len_seq[G->nb_seq] = 0;
	  // get the sequence name ("CHROMOSOME" column in VCF) without the first char (>) nor the last (\n)
	  BUF[strlen(BUF)-1] = '\0';
	  G->name[G->nb_seq] = strndup(&BUF[1], strlen(BUF)-1);
	  G->nb_seq++;
	}
      else
	{
	  k = strlen(BUF);
	  for (i=0; i<k-1; i++)
	    {
	      if (BUF[i]!='N')
		{
		  j = (int) BUF[i];
		  G->data[x++] = (j>>1)&3; // A -> 0, C -> 1, G -> 3, T -> 2
		}
	      else
		{
		  G->data[x++] = 4; // doit etre considere comme zone faible complexite
		}
	      G->len_seq[G->nb_seq-1]++;
	    }
	  }
    }
  fclose(fgen);


  // detection des zones de faible complexite
  if (0)
  {
    int ns;
    long is, lx;
    int v, x, k;
    
    for (ns=0; ns < G->nb_seq; ns++)
      {
	lx = G->pt_seq[ns]; // lx = position dans GENOME du 1er caractere de la sequence numero ns
	for (is=0; is < G->len_seq[ns]-SIZE_SEED; is++)
	  {
	    v=code_seed(&G->data[lx+is]);
	    x = 0;
	    for (k=0; k<16; k++)
	      {
		if (v==x)
		  {
		    G->data[lx+is] = G->data[lx+is]|4;  // marque zone de faible complexite
		    break;
		  }
		x = x + 0x111111;
	      }
	  }
      }
  }

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
