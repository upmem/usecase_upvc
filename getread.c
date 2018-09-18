#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

#define MAX_SIZE_READ 512
#define MAX_SIZE_COM  1024

char COMMENT[MAX_SIZE_COM];
char PEBUF[MAX_SIZE_READ];
int INVNT[4] = {2,3,0,1};

// on retourne un read et son complement 
// la taille est un multiple de 4
// on tronque si c'est un peu plus long
int getSeqFastAQ(FILE *ff, int8_t *read1, int8_t *read2)
{
  int i;
  int offset = 0;
  char c1, c2;
  if (fgets(COMMENT,MAX_SIZE_COM,ff)==NULL) return -1; //lecture commentaire
  if (fgets(PEBUF,MAX_SIZE_READ,ff)==NULL) return -1;  //lecture sequence
  // un commentaire du style >>12 indique qu'il faut omettre les premiers 12 caracteres du read
  if (COMMENT[1] == '>')
    {
      sscanf(COMMENT,"%c%c%d",&c1,&c2,&offset);
    }
  for (i=0; i<SIZE_READ-offset; i++) 
    {
      read1[i] = (((int)PEBUF[i]) >> 1)&3;
      read2[SIZE_READ-i-1-offset] = INVNT[read1[i]];
    }
  for (i=SIZE_READ-offset; i<SIZE_READ; i++) { read1[i] = 0; read2[i] = 0; }
  if (COMMENT[0] == '>') return SIZE_READ;
  if (fgets(COMMENT,MAX_SIZE_COM,ff)==NULL) return -1; //lecture commentaire
  if (fgets(PEBUF,MAX_SIZE_READ,ff)==NULL) return -1; //lecture qualité
  return SIZE_READ;
}


// memorisation dans le buffer BUF_READ des reads paires et de leur complement
// retourne le nombre de reads lu * 2 
int getPEreads(FILE *fpe1, FILE *fpe2, int8_t *BUF_READ, TIMES *CT)
{
  double t1, t2;
  int nb_read = 0;
  t1 = my_clock();
  while (nb_read < MAX_BUF_READ)
    {
      if ((getSeqFastAQ(fpe1,&BUF_READ[nb_read*SIZE_READ],&BUF_READ[(nb_read+1)*SIZE_READ])<=0) 
	  || 
	  (getSeqFastAQ(fpe2,&BUF_READ[(nb_read+2)*SIZE_READ],&BUF_READ[(nb_read+3)*SIZE_READ]) <= 0)
	  ) break;
      nb_read += 4;
    }
  t2 = my_clock();
  CT->get_reads = t2-t1;
  CT->tot_get_reads += t2-t1;
  return nb_read;
}

int getSizeRead(char *name)
{
  FILE *ff;
  int size;
  char filename[1024];

  sprintf(filename,"%s_PE1.fastq",name);
  ff = fopen(filename,"r");
  if (fgets(PEBUF,MAX_SIZE_READ,ff)==NULL) return -1; //lecture commentaire
  if (fgets(PEBUF,MAX_SIZE_READ,ff)==NULL) return -1;  //lecture sequence
  size = strlen(PEBUF) -1;
  size = size / 4;
  size = size * 4;
  fclose(ff);
  return size;
}
