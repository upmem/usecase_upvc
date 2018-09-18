#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"




int printTree(VARTREE *root, int8_t *MAPCOV, FILE *fvcf, FILE *fvar)
{
  VARIANT *v;
  int nb_var = 0;

  if (root == NULL)
    return 0;

  nb_var = printTree(root->left, MAPCOV, fvcf, fvar);
  if (root->left != NULL)
    free(root->left);

  v = root->vars;
  if (v->depth>6)
    {
      fprintf(fvcf, "%s\t%d\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d\n", v->chr, v->pos+1-v->offset, v->ref, v->alt, v->depth,MAPCOV[v->pos]);
      nb_var++;
    }

  if (v->depth>2)
    {
      fprintf(fvar, "%s\t%d\t%s\t%s\ttDEPTH=%d\tCOV=%d\n", v->chr, v->pos+1-v->offset, v->ref, v->alt, v->depth,MAPCOV[v->pos]);
    }

  free(root->vars);

  nb_var += printTree(root->right, MAPCOV, fvcf, fvar);
  if (root->right != NULL)
    free(root->right);

  return nb_var;
}

int printVariants(LIST_VARIANTS *LV, int8_t *MAPCOV, FILE *fvcf, FILE *fvar)
{
  int nb_var;
  nb_var = printTree(LV->vt, MAPCOV, fvcf, fvar);
  return nb_var;
}



int create_vcf(char *name_chr, GENOME *RG, LIST_VARIANTS *LV, int *SUB, int8_t *MAPCOV, TIMES *CT)
{
  double t1, t2;
  FILE *fvcf;
  FILE *fvar;
  char NT[4] = {'A','C','T','G'};
  int ns, is, startPos, i, v, x, cnt, nb_var;
  VARIANT *newvar;
  char filename[1024];

  t1 = my_clock();

  sprintf(filename, "%svars_upvc.vcf", name_chr);
  fvcf = fopen(filename,"w");

  sprintf(filename, "%svars_upvc.var", name_chr);
  fvar = fopen(filename,"w");

  // ####### START OF HEADER #######

  // print vcf version (required)
  fprintf(fvcf, "##fileformat=VCFv4.3\n");

  // print source of VCF file (this program)
  fprintf(fvcf, "##source=UPVC %s\n", VERSION);

  // get the file date
  char filedate[10];
  time_t mytime = time(NULL);
  strftime(filedate, 100, "%Y%d%m", localtime(&mytime));
  fprintf(fvcf, "##fileDate=%s\n", filedate);

  // print reference genome file name
  fprintf(fvcf, "##reference=%s\n", RG->filename);

  // print the column names (fields are tab-delimited in VCF)
  fprintf(fvcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");


  // ####### END OF HEADER #######

  // insert substitution variant in the VARIANT tree
  // for each sequence in the genome
  for (ns=0; ns < RG->nb_seq; ns++)
  {
      startPos = RG->pt_seq[ns];
      // for each position in the sequence
      for (is=0; is<RG->len_seq[ns]; is++)
      {
	// get the substitution at the `startPos+is`th position in the genome
	x = startPos+is;
	v = SUB[x];
	for (i=0; i<4; i++)
	  {
	    // get number of substitutions for the base nt[i] at that position
	    cnt = (v>>(i*8))&0xFF;
	    if (cnt > 7)
	      {
		newvar = (VARIANT*) malloc(sizeof(VARIANT));
		newvar->offset = startPos;
		newvar->chr = RG->name[ns]; 
		newvar->depth = cnt;
		newvar->pos = x;
		newvar->ref[0] = NT[RG->data[x]&3];
		newvar->ref[1] = '\0';
		newvar->alt[0] = NT[i];
		newvar->alt[1] = '\0';
		insertVariants(LV, newvar);
	      }
	  }
      }
  }

  //visuVariants(LV); printf ("\n");

  nb_var = printVariants(LV,MAPCOV,fvcf,fvar);

  fclose(fvcf);
  t2 = my_clock();
  CT->vcf = t2-t1;
  return nb_var;
}

