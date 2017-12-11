#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "upvc.h"

int min (int a, int b) { if (a<b) return a; else return b; }


void displayDPDMatrix(int M[SIZE_NBR4+1][SIZE_NBR4+1], int start)
{
  int i, j;

  for (i=0; i<SIZE_NBR4/2; i++)
    {
      for (j=0; j<SIZE_NBR4/2; j++)
	{
	    printf ("%3d ",M[i][j]);
	}
      printf ("\n");
    }
  printf ("---\n");
  for (i=SIZE_NBR4/2; i<SIZE_NBR4+1; i++)
    {
      for (j=SIZE_NBR4/2; j<SIZE_NBR4+1; j++)
	{
	    printf ("%3d ",M[i][j]);
	}
      printf ("\n");
    }
}

// calcul d'une distance d'alignement par programmation dynamique
// fonction utilisee pour les tests uniquement

void DP(int8_t *s1, int8_t *s2)
{
  int D[SIZE_NBR4][SIZE_NBR4];				
  int P[SIZE_NBR4][SIZE_NBR4];				
  int Q[SIZE_NBR4][SIZE_NBR4];				
  int i, j, d,  QP;

  for (i=0; i<SIZE_NBR4; i++) { D[i][0] = i*COST_SUB; P[i][0] = 99;  Q[i][0] = 99; }
  for (j=0; j<SIZE_NBR4; j++) { D[0][j] = j*COST_SUB; P[0][j] = 99;  Q[0][j] = 99;}

  for (i=1; i<SIZE_NBR4; i++)
    {
      for (j=1; j<SIZE_NBR4; j++)
	{
	  P[i][j] = min(D[i][j-1]+COST_GAPO,P[i][j-1]+COST_GAPE);
	  Q[i][j] = min(D[i-1][j]+COST_GAPO,Q[i-1][j]+COST_GAPE);
	  QP = min(P[i][j],Q[i][j]);
	  d = D[i-1][j-1];
	  if (s1[i-1]!=s2[j-1]) d+=COST_SUB;
	  D[i][j] = min(d,QP);
	}
    }
}


// calcul d'une distance d'alignement par programmation dynamique
// sur les diagonales de la matrice
// on retourne le chemin de backtrack

int DPD(int8_t *s1, int8_t *s2, BACKTRACK *bt)
{
  int D[SIZE_NBR4+1][SIZE_NBR4+1];				
  int P[SIZE_NBR4+1][SIZE_NBR4+1];				
  int Q[SIZE_NBR4+1][SIZE_NBR4+1];				
  int T[SIZE_NBR4+1][SIZE_NBR4+1];				

  int i, ii, j, jj, d,  k, QP, min_score;

  jj = 0;
  ii = 0;

  // pas utile : uniquement pour l'affichage des resultats
  for (i=0; i<SIZE_NBR4+1; i++) for (j=0; j<SIZE_NBR4+1; j++) D[i][j] = 0;

  for (i=0; i<=NB_DIAG/2+1; i++) { P[i][0]=99; Q[i][0]=99; P[0][i]=99; Q[0][i]=99; }
  for (i=0; i<=NB_DIAG/2+1; i++) { D[i][0] = i*COST_SUB; }
  for (j=0; j<=NB_DIAG/2+1; j++) { D[0][j] = j*COST_SUB; }

  for (i=1; i<NB_DIAG/2+1; i++)
    {
      for (j=1; j<=i+NB_DIAG/2; j++)
	{
	  if (D[i][j-1]+COST_GAPO<P[i][j-1]+COST_GAPE) { P[i][j] = D[i][j-1]+COST_GAPO; } else { P[i][j] = P[i][j-1]+COST_GAPE; }
	  if (D[i-1][j]+COST_GAPO<Q[i-1][j]+COST_GAPE) { Q[i][j] = D[i-1][j]+COST_GAPO; } else { Q[i][j] = Q[i-1][j]+COST_GAPE; }
	  if (P[i][j] < Q[i][j]) { QP =  P[i][j]; T[i][j] = 1; } else { QP = Q[i][j]; T[i][j] = 2; }
	  d = D[i-1][j-1];
	  if (s1[i-1]!=s2[j-1]) d+=COST_SUB;
	  if (d<QP) { D[i][j] =  d; T[i][j] = 0; } else { D[i][j] =  QP; }
	}
      Q[i][j]=99; D[i][j]=99;
    }
  for (i=NB_DIAG/2+1; i<SIZE_NBR4-NB_DIAG/2; i++)
    {
      j=i-NB_DIAG/2-1;
      P[i][j] = 99; D[i][j] = 99;
      for (j=i-NB_DIAG/2; j<=i+NB_DIAG/2; j++)
	{
	  if (D[i][j-1]+COST_GAPO<P[i][j-1]+COST_GAPE) { P[i][j] = D[i][j-1]+COST_GAPO; } else { P[i][j] = P[i][j-1]+COST_GAPE; }
	  if (D[i-1][j]+COST_GAPO<Q[i-1][j]+COST_GAPE) { Q[i][j] = D[i-1][j]+COST_GAPO; } else { Q[i][j] = Q[i-1][j]+COST_GAPE; }
	  if (P[i][j] < Q[i][j]) { QP =  P[i][j]; T[i][j] = 1; } else { QP = Q[i][j]; T[i][j] = 2; }
	  d = D[i-1][j-1];
	  if (s1[i-1]!=s2[j-1]) d+=COST_SUB;
	  if (d<QP) { D[i][j] =  d; T[i][j] = 0; } else { D[i][j] =  QP; }
	}
      Q[i][j]=99; D[i][j]=99;
    }
  min_score = 99;
  for (i=SIZE_NBR4-NB_DIAG/2; i<SIZE_NBR4+1; i++)
    {
      j=i-NB_DIAG/2-1;
      P[i][j] = 99; D[i][j] = 99;
      for (j=i-NB_DIAG/2; j<SIZE_NBR4+1; j++)
	{
	  if (D[i][j-1]+COST_GAPO<P[i][j-1]+COST_GAPE) { P[i][j] = D[i][j-1]+COST_GAPO; } else { P[i][j] = P[i][j-1]+COST_GAPE; }
	  if (D[i-1][j]+COST_GAPO<Q[i-1][j]+COST_GAPE) { Q[i][j] = D[i-1][j]+COST_GAPO; } else { Q[i][j] = Q[i-1][j]+COST_GAPE; }
	  if (P[i][j] < Q[i][j]) { QP =  P[i][j]; T[i][j] = 1; } else { QP = Q[i][j]; T[i][j] = 2; }
	  d = D[i-1][j-1];
	  if (s1[i-1]!=s2[j-1]) d+=COST_SUB;
	  if (d<QP) { D[i][j] =  d; T[i][j] = 0; } else { D[i][j] =  QP; }
	}
      if (D[i][SIZE_NBR4]<min_score) { min_score = D[i][SIZE_NBR4]; ii=i; jj=SIZE_NBR4; }
    }
  i=SIZE_NBR4;
  for (j=i-NB_DIAG/2; j<SIZE_NBR4+1; j++)
    if (D[i][j]<min_score) { min_score = D[i][j]; ii=i; jj=j; }

  // backtracking
  // dans bt on met les positions singulieres du chemin
  i = ii;
  j = jj;
  k = 0; 

  // suppression des indels en fin d'alignements
  // T[i][j] contient la trace du chemin
  k = 1;
  while (k==1)
    {
      if (T[i][j] != 0)
	{
	  if (T[i][j] == 1)
	    {
	      j--;
	    }
	  else
	    {
	      i--;
	    }	  
	}
      else
	{
	  k = 0;
	}
    }

  // les conditions i>1 et j>1 suppriment les insertions/deletions
  // en debut d'alignements.
  bt[0].type = CODE_END;
  k = 1;
  while ((i>1) && (j>1))
    {
      if (T[i][j] == 0)
	{
	  i--; j--;
	  if (D[i][j] != D[i-1][j-1]) // substitution
	    {
	      bt[k].type = CODE_SUB;
	      bt[k].ix = i;
	      bt[k].jx = j;  
	      k = k+1;
	    }
	}
      else
	{
	  if (T[i][j] == 1) 
	    {
	      j--;
	      bt[k].type = CODE_INS;
	      bt[k].ix = i;
	      bt[k].jx = j;  
	      k = k+1;
	    }
	  else 
	    {
	      i--;
	      bt[k].type = CODE_DEL;
	      bt[k].ix = i;
	      bt[k].jx = j;  
	      k = k+1;
	    }
	}
    }
  return k;
}

