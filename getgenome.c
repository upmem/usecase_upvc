#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>

#include "upvc.h"


GENOME *get_genome(char* name, TIMES *CT)
{
  
  FILE *fgen;          // fichier genome
  long sizefile;       // taille du fichier
  int  i;
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
	  size_t k = strlen(BUF);
	  for (i=0; i<k-1; i++)
	    {
          int c = toupper(BUF[i]);
	      if (c!='N')
		{
		  G->data[x++] = (c>>1)&3; // A -> 0, C -> 1, G -> 3, T -> 2
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
  {
    int i,ns;
    long is, lx;
    int C[16];
    int8_t k1, k2;
    int x, k;

    int v=0;
    
    for (ns=0; ns < G->nb_seq; ns++)
      {
	lx = G->pt_seq[ns]; // lx = position dans GENOME du 1er caractere de la sequence numero ns
	for (i=0; i<16; i++) C[i] = 0;
	for (is=0; is < SIZE_WINDOW_CPLX; is++)
	  {
	    k1 = (G->data[lx+is]&3) + ((G->data[lx+is+1]&3)<<2);
	    C[k1]++;
	    k2 = (G->data[lx+is+1]&3) + ((G->data[lx+is]&3)<<2);
	    if (k1!=k2) C[k2]++;
	  }
	for (is=0; is < G->len_seq[ns]-SIZE_SEED; is++)
	  {
	    v=code_seed(&G->data[lx+is]); // v contient le code de la graine
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
	    continue;


	    k1 = (G->data[lx+is]&3) + ((G->data[lx+is+1]&3)<<2);
	    C[k1]--;
	    k2 = (G->data[lx+is+1]&3) + ((G->data[lx+is]&3)<<2);
	    if (k1!=k2) C[k2]--;
	    k1 = (G->data[lx+is+SIZE_WINDOW_CPLX]&3) + ((G->data[lx+is+SIZE_WINDOW_CPLX+1]&3)<<2);
	    C[k1]++;
	    k2 = (G->data[lx+is+SIZE_WINDOW_CPLX+1]&3) + ((G->data[lx+is+SIZE_WINDOW_CPLX]&3)<<2);
	    if (k1!=k2)C[k2]++;
	    if (C[k1] >= SIZE_WINDOW_CPLX)
	      {
		if (G->data[lx+is]<4) v++;
		G->data[lx+is] = G->data[lx+is]|4;  // marque zone de faible complexite
	        //for (v=0; v<SIZE_WINDOW_CPLX; v++) printf ("%d",(int)G->data[lx+is+v]); printf (" %d %ld\n",C[k],lx+is);
	      }
	  }
      }
    printf ("v=%d\n",v);
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
