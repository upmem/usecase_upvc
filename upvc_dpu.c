#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

// organisation de la memoire d'un DPU
typedef struct {
                              // index
  int8_t *neighbor_idx;       // tableau qui contient le voisinage (code) des graines du genome de reference
  long *coordinate;           // tableau qui contient le numeros de sequence et l'indice dans a sequence de la graine

                              // pool de reads
  int8_t *neighbor_read;      // tableau qui contient voisinage (code) du read
  int *offset;                // indique pour chaque read l'adresse du 1er voisinage
  int *count;                 // indique pour chaque read le nombre de voisinages
  int *num;                   // indique le numero des reads

                              // resultat
  int *out_score;             // resultat : score
  int *out_num;               // resultat : no du read ou on a trouve un match
  long *out_coord;            // resultat : coordonnees sur le genome ou le read a matche
} MEM_DPU;

MEM_DPU *MDPU; 

// pour le calcul des statistiques

long stat_nb_read[NB_DPU];   // data for ploting the number of reads dispatched by DPU
long stat_nb_nbr[NB_DPU];    // data for ploting the number of distances computed by DPU



int mini (int a, int b) { if (a<b) return a; else return b; }

// calcul d'une distance d'alignement par programmation dynamique
// sur les diagonales de la matrice
// optimisation des resources 
// s'arrete quand le score a depasse le seuil max_score

int ODPD(int8_t *s1, int8_t *s2, int max_score)
{
  int D[2][SIZE_NBR4+1];				
  int P[2][SIZE_NBR4+1];				
  int Q[2][SIZE_NBR4+1];				
  int i, j, d,  lp, pp, QP, min_score;

  for (j=0; j<=NB_DIAG/2+1; j++) { P[0][j]=99; Q[0][j]=99; } P[1][0]=99; Q[1][0]=99; 
  for (j=0; j<=NB_DIAG/2+1; j++) { D[0][j] = j*COST_SUB; } 

  for (i=1; i<NB_DIAG/2+1; i++)
    {
      min_score = 99; // ajout v1.2
      pp = i%2;
      lp = (i-1)%2;
      D[pp][0] = i*COST_SUB;
      for (j=1; j<=i+NB_DIAG/2; j++)
	{
	  P[pp][j] = mini(D[pp][j-1]+COST_GAPO,P[pp][j-1]+COST_GAPE);
	  Q[pp][j] = mini(D[lp][j]+COST_GAPO,Q[lp][j]+COST_GAPE);
	  QP = mini(P[pp][j],Q[pp][j]);
	  d = D[lp][j-1];
	  //if (s1[i-1]!=s2[j-1]) d+=COST_SUB; // delete v1.2
	  if (  (( s1[(i-1)/4] >> (2*((i-1)%4)) )&3) != (( s2[(j-1)/4] >> (2*((j-1)%4)) )&3) ) d+=COST_SUB;
	  D[pp][j] = mini(d,QP);
	  if (D[pp][j] < min_score) min_score = D[pp][j]; // ajout v1.2
	}
      Q[pp][j]=99; D[pp][j]=99;
      if (min_score > max_score) return min_score; // ajout v1.2
    }

  for (i=NB_DIAG/2+1; i<SIZE_NBR4-NB_DIAG/2; i++)
    {
      min_score = 99;
      pp = i%2;
      lp = (i-1)%2;
      j=i-NB_DIAG/2-1;
      P[pp][j] = 99; D[pp][j] = 99;
      for (j=i-NB_DIAG/2; j<=i+NB_DIAG/2; j++)
	{
	  P[pp][j] = mini(D[pp][j-1]+COST_GAPO,P[pp][j-1]+COST_GAPE);
	  Q[pp][j] = mini(D[lp][j]+COST_GAPO,Q[lp][j]+COST_GAPE);
	  QP = mini(P[pp][j],Q[pp][j]);
	  d = D[lp][j-1];
	  //if (s1[i-1]!=s2[j-1]) d+=COST_SUB; // delete v1.2
	  if (   (( s1[(i-1)/4] >> (2*((i-1)%4)) )&3) != (( s2[(j-1)/4] >> (2*((j-1)%4)) )&3) ) d+=COST_SUB;
	  D[pp][j] = mini(d,QP);
	  if (D[pp][j] < min_score) min_score = D[pp][j];
	}
      Q[pp][j]=99; D[pp][j]=99;
      if (min_score > max_score) return min_score;
    }
  min_score = 99;
  for (i=SIZE_NBR4-NB_DIAG/2; i<SIZE_NBR4+1; i++)
    {
      pp = i%2;
      lp = (i-1)%2;
      j=i-NB_DIAG/2-1;
      P[pp][j] = 99; D[pp][j] = 99;
      for (j=i-NB_DIAG/2; j<SIZE_NBR4+1; j++)
	{
	  P[pp][j] = mini(D[pp][j-1]+COST_GAPO,P[pp][j-1]+COST_GAPE);
	  Q[pp][j] = mini(D[lp][j]+COST_GAPO,Q[lp][j]+COST_GAPE);
	  QP = mini(P[pp][j],Q[pp][j]);
	  d = D[lp][j-1];
	  //if (s1[i-1]!=s2[j-1]) d+=COST_SUB; // delete v1.2
	  if (   (( s1[(i-1)/4] >> (2*((i-1)%4)) )&3) != (( s2[(j-1)/4] >> (2*((j-1)%4)) )&3) ) d+=COST_SUB;
	  D[pp][j] = mini(d,QP);
	}
      if (D[pp][SIZE_NBR4]<min_score) min_score = D[pp][SIZE_NBR4];
    }
  i=SIZE_NBR4;
  pp = i%2;
  for (j=i-NB_DIAG/2; j<SIZE_NBR4+1; j++)
    if (D[pp][j]<min_score) min_score = D[pp][j];
  return min_score;
}

int TT[256] = {   0 , 10 , 10 , 10 , 10 , 20 , 20 , 20 , 10 , 20 , 20 , 20 , 10 , 20 , 20 , 20 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 10 , 20 , 20 , 20 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 , 20 , 30 , 30 , 30 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 ,
		 20 , 30 , 30 , 30 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 , 30 , 40 , 40 , 40 };


// version optimisee de noDP
// retourne un score si pas de détection d'indels
// sinon (détection d'indels) retourne -1. Dans ce cas on lancera la procédure ODPD
// fonction qui utilise le casting de pointeurs

int noDP(int8_t *s1, int8_t *s2, int max_score)
{
  int i, j, x, v, V, V1, V2;
  int *X1, *X2;
  int score = 0;
  for (i=0; i<SIZE_NBR; i++)
    {
      x = (int) (s1[i]^s2[i]);
      x = x&0xFF;
      v = TT[x];
      if (v > COST_SUB) // indique si on a plus d'une difference
	{
	  j=i+1;
	  if (j<SIZE_NBR-3)
	    {
	      X1 =  (int *) (&s1[j]);
	      X2 =  (int *) (&s2[j]);
	      V1 = *X1;
	      V2 = *X2;
	      V = (V1 ^ V2)&0xFFFFFF; // on regarde si les 8 caracteres suivants sont identiques
	      if (V!=0) // si caracteres differents on test les indels
		{
		  V = (V1 ^ (V2>>2))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V1 ^ (V2>>4))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V1 ^ (V2>>6))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V1 ^ (V2>>8))&0xFFFFFF;
		  if (V==0) return -1;
		  
		  V = (V2 ^ (V1>>2))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V2 ^ (V1>>4))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V2 ^ (V1>>6))&0xFFFFFF;
		  if (V==0) return -1;
		  V = (V2 ^ (V1>>8))&0xFFFFFF;
		  if (V==0) return -1;
		}
	    }
	}
      score += v;
      if (score > max_score) break;
    }
  return score;
}

void print2NBR(int8_t *s1, int8_t *s2)
{
  int i,j;
  int8_t x, x1, x2;
  for (i=0; i<SIZE_NBR; i++)
    {
      for (j=0; j<4; j++)
	{
	  x = (s1[i]>>(2*j))&3;
	  printf ("%x",x);
	}
    }
  printf ("\n");
  for (i=0; i<SIZE_NBR; i++)
    {
      for (j=0; j<4; j++)
	{
	  x1 = (s1[i]>>(2*j))&3;
	  x2 = (s2[i]>>(2*j))&3;
	  if (x1!=x2) printf ("|"); else printf (" ");
	}
    }
  printf ("\n");
  for (i=0; i<SIZE_NBR; i++)
    {
      for (j=0; j<4; j++)
	{
	  x = (s2[i]>>(2*j))&3;
	  printf ("%x",x);
	}
    }
  printf ("\n\n");
}





// mapping des reads
// M point sur une memoire d'un des DPUs
// M.num[nr] == -1 : ==> marqueur de fin

void align(int numdpu)
{
  int nr, ix, score_noDP, score_ODPD, score, offset, mini, nb_map_start;
  int nb_map = 0;
  MEM_DPU M = MDPU[numdpu];

  stat_nb_read[numdpu]=0;
  stat_nb_nbr[numdpu]=0;

  nr = 0;
  M.out_num[nb_map] = -1;
  while (M.num[nr] != -1)   // nr = indice qui parcours les reads a traiter
    {
      //stat_nb_read[numdpu]++;
      offset = M.offset[nr];                            // offset = adresse du 1er voisinage
      mini = MAX_SCORE;
      nb_map_start = nb_map;
      for (ix=0; ix<M.count[nr]; ix++)                  // count = nombre de voisinages
	{
	  // stat_nb_nbr[numdpu]++;
	  score_noDP = noDP(&M.neighbor_read[nr*SIZE_NBR],&M.neighbor_idx[(offset+ix)*SIZE_NBR],mini);
	  score = score_noDP;
	  if (score_noDP == -1)
	    {
	      score_ODPD = ODPD(&M.neighbor_read[nr*SIZE_NBR],&M.neighbor_idx[(offset+ix)*SIZE_NBR],mini);
	      score = score_ODPD;
	    }
	  if (score <= mini)
	    {
	      if (score <  mini)
		{
		  mini = score;
		  nb_map = nb_map_start;
		}
	      if (nb_map < MAX_ALIGN-1)
		{
		  M.out_num[nb_map] = M.num[nr];
		  M.out_coord[nb_map] = M.coordinate[offset+ix];
		  M.out_score[nb_map] = score;
		  nb_map++;
		  M.out_num[nb_map] = -1;
		}
	    }
	}
      nr++;
    }
  //printf ("STAT %d %ld %ld\n",numdpu,stat_nb_read[numdpu],stat_nb_nbr[numdpu]);
}

// allocation globale de la memoire DPU
void malloc_dpu()
{
  int d;
  MDPU = (MEM_DPU *) malloc(sizeof(MEM_DPU)*NB_DPU);
  for (d=0; d<NB_DPU; d++)
    {
      MDPU[d].neighbor_read  = (int8_t *) malloc(sizeof(int8_t)*SIZE_NBR*MAX_NB_DPU_READ); 
      MDPU[d].count          = (int *)    malloc(sizeof(int)*MAX_NB_DPU_READ); 
      MDPU[d].offset         = (int *)    malloc(sizeof(int)*MAX_NB_DPU_READ); 
      MDPU[d].num            = (int *)    malloc(sizeof(int)*MAX_NB_DPU_READ); 
      MDPU[d].out_num        = (int *)    malloc(sizeof(int)*MAX_ALIGN); 
      MDPU[d].out_coord      = (long *)   malloc(sizeof(long)*MAX_ALIGN); 
      MDPU[d].out_score      = (int  *)   malloc(sizeof(long)*MAX_ALIGN); 
    }
}

// allocation de la memoire des index
void malloc_neighbor_idx (int d, int k) 
{
  MDPU[d].neighbor_idx  = (int8_t *) malloc(sizeof(int8_t)*k*SIZE_NBR); 
  MDPU[d].coordinate    = (long *) malloc(sizeof(long)*k); 
}

// liberation de la memoire
void free_dpu()
{
  int i;
  for (i=0; i<NB_DPU; i++)
    {
      free (MDPU[i].neighbor_idx);
      free (MDPU[i].neighbor_read);
      free (MDPU[i].coordinate);
      free (MDPU[i].count);
      free (MDPU[i].offset);
      free (MDPU[i].num);
      free (MDPU[i].out_num);
      free (MDPU[i].out_score);
      free (MDPU[i].out_coord);
    }
  free(MDPU);
}

void write_neighbor_read (int d, int k, int8_t *v) 
{ 
  int i;
  for (i=0; i<SIZE_NBR; i++)
    MDPU[d].neighbor_read[k*SIZE_NBR+i] = v[i];
}

void write_neighbor_idx (int d, int k, int8_t *v) 
{ 
  int i;
  for (i=0; i<SIZE_NBR; i++)
    MDPU[d].neighbor_idx[k*SIZE_NBR+i] = v[i];
}

void write_coordinate  (int d, int k, long v) { MDPU[d].coordinate[k] = v; }
void write_count       (int d, int k, int v) { MDPU[d].count[k] = v; }
void write_offset      (int d, int k, int v) { MDPU[d].offset[k] = v; }
void write_num         (int d, int k, int v) { MDPU[d].num[k] = v; }

int read_out_num (int d, int k) { return MDPU[d].out_num[k]; }
int read_out_score (int d, int k) { return MDPU[d].out_score[k]; }
long read_out_coord (int d, int k) { return MDPU[d].out_coord[k]; }



