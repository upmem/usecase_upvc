#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

int cmpalign(const void * a, const void * b) 
{
  ALIGN *A = (ALIGN *) a;
  ALIGN *B = (ALIGN *) b;
  if (A->num_read == B->num_read) 
    {
      if (A->score > B->score)
	return 1;
      else
	return -1;
    }
  if (A->num_read > B->num_read)
    return 1;
  else
    return -1;
}


void getSizeInsert (ALIGN * A, int nb_match) 
{
  SIZE_INSERT_MEAN = 400;
  SIZE_INSERT_STD = 50;
}

int check_pair(ALIGN A1, ALIGN A2) 
{
  int pos1, pos2;
  int ret = 0;
  int size_insert_min = SIZE_INSERT_MEAN - 3*SIZE_INSERT_STD;
  int size_insert_max = SIZE_INSERT_MEAN + 3*SIZE_INSERT_STD;
  
  if (A1.coord_seq == A2.coord_seq) 
    { // make sure the two reads align to the same reference sequence
      pos1 = A1.coord_pos;
      pos2 = A2.coord_pos;

      // verification de l'ecartement entre reads
      if ((abs(pos2 - pos1 + SIZE_READ) > size_insert_min) && (abs(pos2 - pos1 + SIZE_READ) < size_insert_max)) 
	{
	  ret = 1;
	}
    }
  return ret;
}

// codage des alignements
// le code est recupere dans la chaine de caracteres codali
// c'est une suite d'alignements elementaires 
// codage d'une substitution : S pos x   [ pos = position dans le read de la substitution, x = A,C,G ou T ]
// codage d'une insertion :    I pos x+  [ pos = position dans le read de l'insertion, x+ =  1 ou plusieurs caracteres A,C,G ou T ]
// codage d'une deletion :     D pos x+  [ pos = position dans le read de la deletion, x+ =  1 ou plusieurs caracteres A,C,G ou T ]
// le dernier caratere est le symbole X
// example S 12 A D 56 A T G I 87 T C X ==> substitution (A) position 12, deletion (ATG) position 56, insertion (TC) position 87
// retourne la taille de la chaine de caracteres

int code_alignment(int8_t *code, int score, int8_t *gen, int8_t *read) 
{
  int i, k, s, x, ix, jx;
  BACKTRACK BT[SIZE_READ];

  // codage substitution          CODE_SUB pos x
  //        deletion              CODE_DEL pos x+
  //        insertion             CODE_INS pos x+
  //        fin                   CODE_END
  //
  // x  = A | C | G | T
  // x+ = une sequence d'au moins 1 X (i.e. A, C, G ou T)
  // pos = entier (8 bits) : donne la position du variant par rapport au debut du read

  if (score == 0) 
    {
      code[0] = CODE_END;
      return 1;
    }

  // on regarde d'abord si on a a faire a des erreurs de substitution uniquement (c'est vrai la plupart du temps)
  // dans ce cas, on verifie simplement la diagonale de la matrice et on repere les positions de subsitutions
  k = 0;
  s = 0;
  for (i = SIZE_SEED; i < SIZE_NBR4 + SIZE_SEED; i++) 
    {
      if ((gen[i] & 3) != read[i]) 
	{
	  s = s + COST_SUB;
	  code[k++] = CODE_SUB;
	  code[k++] = i;
	  code[k++] = read[i];
	  if (s > score) break;
	}
    }
  code[k++] = CODE_END;
  if (s == score) return k;
		
  // sinon, on recalcule la matrice (sur quelques diagonales uniquement)
  // dans BT on recupere le chemin du backtracking
  // x = taille du chemin

  x = DPD(&gen[SIZE_SEED], &read[SIZE_SEED], BT);

  x--; // indice sur le chemin
  k = 0; // indice sur le code de l'alignement
  while (x > 0) 
    {
      if (BT[x].type == CODE_SUB) // substitution
	{
	  code[k++] = CODE_SUB;
	  code[k++] = BT[x].jx + SIZE_SEED - 1;
	  code[k++] = read[BT[x].jx + SIZE_SEED - 1];
	  x--;
	} 
      else 
	{
	  if (BT[x].type == CODE_DEL) // deletion
	    {
	      code[k++] = CODE_DEL;
	      code[k++] = BT[x].ix + SIZE_SEED;
	      code[k++] = gen[BT[x].ix + SIZE_SEED] & 3;
	      jx = BT[x].jx;
	      x--;
	      while ((BT[x].type == CODE_DEL) && (jx == BT[x].jx)) 
		{
		  code[k++] = gen[BT[x].ix + SIZE_SEED] & 3;
		  x--;
		}
	    } 
	  else 
	    {
	      code[k++] = CODE_INS; // insertion
	      code[k++] = BT[x].jx + SIZE_SEED - 1;
	      code[k++] = read[BT[x].jx + SIZE_SEED];
	      ix = BT[x].ix;
	      x--;
	      while ((BT[x].type == CODE_INS) && (ix == BT[x].ix)) 
		{
		  code[k++] = read[BT[x].jx + SIZE_SEED];
		  x--;
		}
	    }
	}
    }
  code[k++] = CODE_END;
  return k;
}

void set_variant(ALIGN A, GENOME *RG, int8_t *BUF_READ, LIST_VARIANTS *LV, int *SUB, int8_t *MAPCOV) 
{
  int i, k,  type,  snp, ptref, ptalt, lref, lalt, x;
  int8_t codali[256];
  int8_t *READ;
  VARIANT *newvar = (VARIANT *) malloc(sizeof(VARIANT));
  char NT[4] = {'A','C','T','G'};    
  int pos = RG->pt_seq[A.coord_seq] + A.coord_pos; // position dans le genome qui correspond au premier charactere du read mappe

  // mise a jour de MAPCOV = nombre de reads qui mapped a cette position dans le genome
  // a affiner avec les rounds
  for (i=0; i<SIZE_READ; i++) MAPCOV[pos+i] += 1;

  // codage de l'alignement
  // on recupere dans le tableau codali le code des alignements (voir la fonction code_alignment)
  READ = &BUF_READ[A.num_read * SIZE_READ];
  code_alignment(codali, A.score, &RG->data[pos], READ);
			
  // mise a jour du champ cntSNP (compteur de SNPs)
  // et gestion des indels
  i = 0;
  while (codali[i] != CODE_END) 
    {
      type = codali[i];
      ptalt = codali[i+1];         // ptalt pointe sur la position du variant dans le read
      ptref = pos + ptalt;         // ptref pointe sur la position du variant dans le genome
      if (type == CODE_SUB)        // variant de type substitution : S pos x  { codali[i], codali[i+1], codali[i+2] }
	{
	  snp = codali[i+2];       // SNP = 0,1,2,3  (code A,C,T,G)
	  x = 1 << (snp * 8);      // positionne +1 sur le bon compteur
	  SUB[ptref] += x;         // SUB[ptref] = entier 32 bits qui contient 4 compteurs 8 bits
	  i = i+3;
	} 
      else if (type == CODE_INS)    // variant de type insertion : I pos nt+
	{
	  newvar = (VARIANT*) malloc(sizeof(VARIANT));
	  newvar->offset = RG->pt_seq[A.coord_seq];
	  newvar->chr = RG->name[A.coord_seq]; 
	  newvar->depth = 1;
	  lref = ptref;
	  lalt = ptalt;
	  i = i+2;
	  while (codali[i] < 4) // lecture des codes des caracteres A, C, G, T
	    {
	      lalt++;
	      i++;
	    }
	  if ((ptref > 26121800) && (ptref <26121900))
	  {
	    int k;
	    printf ("ptref = %d  ptalt=%d\n",ptref,ptalt);
	    for (k=ptref-20; k<ptref; k++) printf ("%c",NT[RG->data[k]]); printf (" %c ",NT[RG->data[ptref]]);
	    for (k=ptalt; k<lalt; k++) printf (" ");
	    for (k=ptref+1; k<ptref+20; k++) printf ("%c",NT[RG->data[k]]); printf ("\n");
	    for (k=ptalt-20; k<ptalt; k++) printf ("%c",NT[READ[k]]); printf (" %c ",NT[READ[ptalt]]);
	    for (k=ptalt+1; k<ptalt+20; k++) printf ("%c",NT[READ[k]]); printf ("\n");
	    for (k=ptalt-20; k<ptalt; k++) printf (" "); printf (" %c ",NT[READ[ptalt]]);
	    for (k=ptalt+1; k<=lalt; k++) printf ("%c",NT[READ[k]]); printf ("\n");
	  }

	  // normalisation de la position de l'insertion
	  while (RG->data[lref] == READ[lalt]) { lref--; lalt--; ptref--; ptalt--; }

	  if ((ptref > 26121800) && (ptref <26121900))
	  {
	    int k;
	    printf ("ptref = %d\n",ptref);
	    for (k=ptref-20; k<ptref; k++) printf ("%c",NT[RG->data[k]]); printf (" %c ",NT[RG->data[ptref]]);
	    for (k=ptalt; k<lalt; k++) printf (" ");
	    for (k=ptref+1; k<ptref+20; k++) printf ("%c",NT[RG->data[k]]); printf ("\n");
	    for (k=ptalt-20; k<ptalt; k++) printf ("%c",NT[READ[k]]); printf (" %c ",NT[READ[ptalt]]);
	    for (k=ptalt+1; k<ptalt+20; k++) printf ("%c",NT[READ[k]]); printf ("\n");
	    for (k=ptalt-20; k<ptalt; k++) printf(" "); printf (" %c ",NT[READ[ptalt]]);
	    for (k=ptalt+1; k<=lalt; k++) printf ("%c",NT[READ[k]]); printf ("\n");
	  }

	  newvar->pos = ptref;
	  newvar->ref[0] = NT[RG->data[ptref]&3];
	  newvar->ref[1] = '\0';
	  k = 0;
	  while (ptalt <= lalt) { newvar->alt[k] = NT[READ[ptalt]&3]; k++; ptalt++; }
	  newvar->alt[k] = '\0';
	  insertVariants(LV, newvar);
	}
      else if (type == CODE_DEL)
	{
	  newvar = (VARIANT*) malloc(sizeof(VARIANT));
	  newvar->offset = RG->pt_seq[A.coord_seq];
	  newvar->chr = RG->name[A.coord_seq]; 
	  newvar->depth = 1;
	  // ptref--; // difference avec insertion !
	  lref = ptref;
	  lalt = ptalt;
	  i = i+2;
	  while (codali[i] < 4) // lecture des codes des caracteres A, C, G, T
	    {
	      lref++;
	      i++;
	    }
	  // normalisation de la position de la deletion
	  while (RG->data[lref] == READ[lalt]) { lref--; lalt--; ptref--; ptalt--; }
	  newvar->pos = ptref;
	  newvar->alt[0] = NT[RG->data[ptref]]; 
	  newvar->alt[1] = '\0';
	  k = 0;
	  while (ptref <= lref) {newvar->ref[k] = NT[RG->data[ptref]&3]; k++; ptref++; }
	  newvar->ref[k] = '\0';
	  insertVariants(LV, newvar);
	}
    }
}

int process_read(GENOME *RG, int8_t *BUF_READ, int nb_read, LIST_VARIANTS *LV, int *SUB, int8_t *MAPCOV, FILE *f1, FILE *f2, int round, TIMES *CT) 
{
  double t1, t2;
  ALIGN *A = (ALIGN *) malloc(sizeof(ALIGN) * MAX_ALIGN * NB_DPU);
  int N, i, j, k, type, numpair, nb_match, numdpu, num_read, pa, pb;
  int8_t *READ;
  char NT[4] = {'A','C','T','G'};
  long coord;
  int nb_read_map = 0;
  int offset[4], score[4], nbread[4];
  int *TMP = (int *) malloc(sizeof(int) * MAX_ALIGN);

  t1 = my_clock();

  // on recupere dans A les alignements en provenance des DPUs
  nb_match = 0;
  for (numdpu = 0; numdpu < NB_DPU; numdpu++) 
    {
      k = 0;
      while ((num_read = read_out_num(numdpu, k)) != -1) 
	{
	  A[nb_match].num_read = num_read;
	  coord = read_out_coord(numdpu, k);
	  A[nb_match].coord_seq = (int) (coord >> 32);
	  A[nb_match].coord_pos = (int) (coord & 0xFFFFFFFF);
	  A[nb_match].score = read_out_score(numdpu, k);
	  k++;
	  nb_match++;
	}
    }

  // tri des alignements par leur numero et leur score
  qsort(A, nb_match, sizeof(ALIGN), cmpalign);

  // If the mean insert size and standard deviation have not been yet calculated, then calculate them
  if (!SIZE_INSERT_MEAN) 
    {
      getSizeInsert(A, nb_match);
    }

  // le numero d'une paire est donnee par num_read / 4 (voir la fonction dispatch)
  // les reads sont types de par leur numero (voir la fonction dispatch)
  // type =  num_read%4 ==  0 ==> read PE1
  //                        1 ==> complement read PE1
  //                        2 == >read PE2
  //                        3 ==> complement read PE2
  // les paires a considerer sont donc les couples [0,3] et [1,2]

  i = 0;
  while (i < nb_match) 
    {
      // offset[x] = indice dans A du 1er read correspondant au type x
      // nbread[x] = nombre de reads du type x
      // score[x] = score min du type x
      for (k = 0; k < 4; k++) 
	{
	  score[k] = 1000;
	  offset[k] = -1;
	  nbread[k] = 0;
	}
      numpair = A[i].num_read / 4;
      while ((i < nb_match) && (numpair == A[i].num_read / 4)) // on recupres tous les reads d'une paire
	{
	  type = A[i].num_read % 4;
	  // normalement, la fonction cmpalign tri par numero de read et par score
	  // ainsi on ne retient que les reads avec les scores min
	  if (A[i].score < score[type]) 
	    {
	      score[type] = A[i].score;
	      offset[type] = i;
	      nbread[type] = 1;
	    } 
	  else 
	    {
	      if (A[i].score == score[type]) nbread[type]++;
	    }
	  i++;
	}
      N = 0;
      // traitement d'une paire [0,3]
      // on regarde la compatibilité de tous les reads de type 0 avec les reads de type 3
      for (pa = offset[0]; pa < offset[0] + nbread[0]; pa++)
	for (pb = offset[3]; pb < offset[3] + nbread[3]; pb++)
	  if (check_pair(A[pa], A[pb]))
	    if (N < MAX_ALIGN - 2) 
	      {
		TMP[N++] = pa;
		TMP[N++] = pb;
	      }
      // traitement d'une paire [1,2]
      // on regarde la compatibilité de tous les reads de type 1 avec les reads de type 2
      for (pa = offset[1]; pa < offset[1] + nbread[1]; pa++)
	for (pb = offset[2]; pb < offset[2] + nbread[2]; pb++)
	  if (check_pair(A[pb], A[pa]))
	    if (N < MAX_ALIGN - 2) 
	      {
		TMP[N++] = pa;
		TMP[N++] = pb;
	      }
      if (N == 2) // mapping unique
	{
	  set_variant(A[TMP[0]], RG, BUF_READ, LV, SUB, MAPCOV);
	  set_variant(A[TMP[1]], RG, BUF_READ, LV, SUB, MAPCOV);
	  nb_read_map += 2;
	}
      else
	{ // les reads non mappes sont mis dans un fichier en idiquant dans le commentaire
          // l'offset a faire sur les premiers caracteres pour le round suivant
	  READ = &BUF_READ[numpair * SIZE_READ];
	  fprintf (f1,">>%d\n",SIZE_SEED*(round+1));
	  for (j=SIZE_SEED;j<SIZE_READ;j++) fprintf(f1,"%c",NT[READ[j]&3]); 
	  for (j=0; j<SIZE_SEED; j++) fprintf (f1,"A"); fprintf(f1,"\n");

	  READ = &BUF_READ[(numpair+2) * SIZE_READ];
	  fprintf (f2,">>%d\n",SIZE_SEED*(round+1));
	  for (j=SIZE_SEED;j<SIZE_READ;j++) fprintf(f2,"%c",NT[READ[j]&3]); 
	  for (j=0; j<SIZE_SEED; j++) fprintf (f2,"A"); fprintf(f2,"\n");
	}
    }


  free(TMP);
  free(A);
  t2 = my_clock();
  CT->process_read = t2 - t1;
  CT->tot_process_read += t2 - t1;
  return nb_read_map;
}

