#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

// calcul du code d'une graine de taille SIZE_SEED
int code_seed(int8_t *SEQ)
{
  int i;
  int v=0;
  for (i=0; i<SIZE_SEED; i=i+1) 
    {
      if (SEQ[i] >= 4) return -1;
      v=(v<<2)+SEQ[i]; 
    }
  return v;
}

void code_neighbor(int8_t *SEQ, int8_t *CODE)
{
  int i,j;
  for (i=0; i<SIZE_NBR; i++)
    {
      j = i*4;
      CODE[i] = ((SEQ[j+3]&3)<<6) + ((SEQ[j+2]&3)<<4) + ((SEQ[j+1]&3)<<2) + (SEQ[j]&3);
    }
}

void decode_neighbor(int8_t *CODE, int8_t *SEQ)
{
  int i;
  for (i=0; i<SIZE_NBR*4; i=i+4)
    {
      SEQ[i]   = (CODE[i/4])&3;
      SEQ[i+1] = (CODE[i/4]>>2)&3;
      SEQ[i+2] = (CODE[i/4]>>4)&3;
      SEQ[i+3] = (CODE[i/4]>>6)&3;
    }
}

