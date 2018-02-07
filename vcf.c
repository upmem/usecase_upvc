#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

int create_vcf(char *filename, GENOME *RG, VARINDEL **INDEL, int *SUB, int8_t *MAPCOV, TIMES *CT)
{
  double t1, t2;
  FILE *fvcf;
  char nt[4] = {'A','C','T','G'};
  int ns, is, pos, i, j, v, cnt, c, type, size, count, content, nb_var;
  VARINDEL *indel;

  t1 = my_clock();
  fvcf = fopen(filename,"w");

  for (ns=0; ns < RG->nb_seq; ns++)
    {
      pos = RG->pt_seq[ns];
      for (is=0; is<RG->len_seq[ns]; is++)
	{
	  v = SUB[pos+is];
	  j = 0;
	  for (i=0; i<4; i++)
	    {
	      cnt = (v>>(i*8))&0xFF;
	      if (cnt > 5)
		{
		  SUB[pos+is] = CODE_SUB + (1<<4) + (cnt<<8) + (i<<16); j++;
		}
	    }
	  if (j==0) SUB[pos+is] = CODE_END;
	}
    }

  for (i=0; i<SIZE_LIST_INDEL; i++)
    {
      indel = INDEL[i];
      while (indel != NULL)
	{
	  if (indel->count > 5)
	    {
	      SUB[indel->pos_seq] = indel->type + (indel->size<<4) + (indel->count<<8) + (indel->content<<16);
	    }
	  indel = (VARINDEL *) indel->next;
	}
    }

  nb_var = 0;
  for (ns=0; ns < RG->nb_seq; ns++)
    {
      pos = RG->pt_seq[ns];
      for (is=0; is<RG->len_seq[ns]; is++)
	{
	  v = SUB[pos+is];
	  type = v&0xF;
	  if (type != CODE_END)
	    {
	      if (type == CODE_SUB) type = 'S'; else if (type == CODE_INS) type = 'I'; else type = 'D';
	      fprintf (fvcf,"%d %d %c ", ns,is,type);
	      size = (v>>4)&0xF;
	      count = (v>>8)&0xFF;
	      content = v>>16;
	      for (j=size-1; j>=0; j--) fprintf(fvcf,"%c",nt[(content>>(2*j))&3]);
	      c = 0;
	      for (j=is-SIZE_READ; j<is; j++) if (j>0) c+= MAPCOV[j];
	      fprintf (fvcf," (%d/%d)\n",count,c);
	      nb_var++;
	    }
	}
    }

  fclose(fvcf);
  t2 = my_clock();
  CT->vcf = t2-t1;
  return nb_var;
}
