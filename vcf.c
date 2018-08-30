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
  int depth, ns, is, startPos, i, j, v, x, cnt, type, size, content, nb_var, pos, posnt;
  char *ref, *alt, *variation, *chr;
  VARINDEL *indel;
  VARIANT *varlist, *newvar;
  VARTREE* vartree = NULL;

  ref = (char*) calloc(SIZE_ALLELE, sizeof(char));
  alt = (char*) calloc(SIZE_ALLELE, sizeof(char));
  variation = (char*) calloc(SIZE_ALLELE, sizeof(char));

  varlist = (VARIANT*) malloc(sizeof(VARIANT));
  varlist -> next = NULL;

  t1 = my_clock();
  fvcf = fopen(filename,"w");

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
		  j = 0;
		  for (i=0; i<4; i++)
		  {
			// get number of substitutions for the base nt[i] at that position
			cnt = (v>>(i*8))&0xFF;
			// if the number of subs exceeds the threshold of 5, encode the substitution
			// 32 bits are re-divided as follows (from least to most-significant bits):
			// type (4 bits)
			// size (4 bits)
			// count (8 bits)
			// content (16 bits)
			if (cnt > 8)
			{
			  SUB[x] = CODE_SUB + (1<<4) + (cnt<<8) + (i<<16);
			  j++;
			}
		  }
		  // if that position did not have any substitutions greater than 5, then mark it as normal (no variant found)
		  if (j==0) SUB[x] = CODE_END;
      }
  }

  // iterate through the indel linked-list
  for (i=0; i<SIZE_LIST_INDEL; i++)
  {
    indel = INDEL[i];
    while (indel != NULL)
	{
      // if the number of indels exceeds the threshold of 7, encode the indels as per the above template
	  if ((indel->count > 7)&&(MAPCOV[indel->pos_seq]>=20))
	  {
	      SUB[indel->pos_seq] = indel->type + (indel->size<<4) + (indel->count<<8) + (indel->content<<16);
	  }
	  indel = (VARINDEL *) indel->next;
	}
  }

  // counter for the number of variants found
  nb_var = 0;

  // iterate through each variant and output information about it to the VCF file
  for (ns=0; ns < RG->nb_seq; ns++)
  {
    // get the chromosome number from the reference genome metadata
    chr = RG->name[ns];
    // get the starting position of the reference genome (first nucleotide compared)
    startPos = RG->pt_seq[ns];
    // iterate through each nucleotide and see if a variant was found for it
    for (is=0; is<RG->len_seq[ns]; is++)
      {
	// get the encoded variant if it was found
	v = SUB[startPos+is];
	// get the last 4 bits (dedicated to be the variant type) of the encoded variant
	type = v&0xF;
	// get the bits [8-16] as depth
	depth = (v>>8)&0xFF;
	// if it is a variant, process accordingly
	if (type != CODE_END)
	  {
	    // get the type of variant, (S)ubstitution, (I)nsertion, or (D)eletion
	    if (type == CODE_SUB)
	      type = 'S';
	    else if (type == CODE_INS)
	      type = 'I';
	    else
	      type = 'D';

	    // get the size of the variant (always 1 for substitutions)
	    size = (v>>4)&0xF;
	    // get the nucleotide contents of the variant
	    content = v>>16;
	    // get the contents of each variant
	    for (j=size-1, i=0; j>=0; j--, i++)
	      {
		variation[i] = nt[(content>>(2*j))&3];
	      }
	    // terminate the string properly
	    variation[i] = '\0';
		  
	    // if it is a substitution
	    if (type == 'S')
	      {
		// for substitutions, the position is unchanged (no frameshift happened)
		pos = startPos+is;
		// get the nucleotide at that position
		posnt = RG->data[pos];
		// if it is in a weak-spot of the reference genome, then do not count it as a variant, and continue the loop
		if (posnt > 3)
		  continue;
		// otherwise, get the reference and alternate (variant) nucleotides
		sprintf(ref, "%c", nt[posnt]);
		sprintf(alt, "%c", variation[0]);
	      }
	    // if it is an insertion
	    else if (type == 'I')
	      {
		// the reference position is unchanged because the nucleotide is added to its 3' side
		pos = startPos+is;
		// get the nucleotide at that position
		posnt = RG->data[pos];
		// do not process it if it is a weak-spot
		if (posnt > 3)
		  continue;
		// output the reference nucleotide
		sprintf(ref, "%c", nt[posnt]);
		// output the reference nucleotide and append the insertion to it
		sprintf(alt, "%c%s", nt[posnt], variation);
	      }
	    // if the variant is not a substitution or an insertion, it is a deletion
	    else
	      {
		// the reference position is changed because the original nucleotide was deleted, so the reference should be the one before it
		pos = startPos + is - 1;
		// get the nucleotide before the deletion (at its 5' side)
		posnt = RG->data[pos];
		// do not process it if it is a weak-spot
		if (posnt > 3)
		  continue;
		// print the reference position and the portion that was deleted
		sprintf(ref, "%c%s", nt[posnt], variation);
		// print the new content at that position (same as ref minus the deletion)
		sprintf(alt, "%c", nt[posnt]);
	      }
	    
	    if (strcmp(ref,alt) != 0) {
	      newvar = (VARIANT*) malloc(sizeof(VARIANT));
	      newvar->chr = strdup(chr);
	      newvar->offset = startPos;
	      newvar->pos = pos;
	      newvar->ref = strdup(ref);
	      newvar->alt = strdup(alt);
	      newvar->cov = MAPCOV[pos];
	      newvar->depth = depth;
	      newvar->next = NULL;
	      normalize(RG, newvar);
	      vartree = insert(vartree, newvar);
	      nb_var++;
	    }
	  }
      }
  }

  inorder(vartree, fvcf);


  free(vartree);
  free(ref);
  free(alt);
  free(variation);

  fclose(fvcf);
  t2 = my_clock();
  CT->vcf = t2-t1;
  return nb_var;
}

