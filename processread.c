#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "upvc.h"

int cmpalign (const void * a, const void * b)
{
  ALIGN *A = (ALIGN *)a;
  ALIGN *B = (ALIGN *)b;
  if (A->num_read == B->num_read)
    {
      if (A->score > B->score) return 1; else return -1;
    }
  if (A->num_read > B->num_read) return 1; else return -1;
}

int check_pair1(ALIGN A1, ALIGN A2)
{
  int pos1, pos2;
  int ret = 0;
  int size_insert_min = SIZE_READ;
  int size_insert_max = SIZE_READ+200;
  if (A1.coord_seq == A2.coord_seq) // on s'assure que les 2 reads sont sur la meme sequence
    {
      pos1 = A1.coord_pos;
      pos2 = A2.coord_pos;
      if (pos1 > pos2) // verification des orientations des reads paires
	{
	  if (((pos1 - pos2) > size_insert_min) && ((pos1 - pos2) < size_insert_max)) // verification de l'ecartement entre reads
	    {
	      ret = 1;
	    }
	}
    }
  return ret;
}

int check_pair2(ALIGN A1, ALIGN A2)
{
  int pos1, pos2;
  int ret = 0;
  int size_insert_min = SIZE_READ;
  int size_insert_max = SIZE_READ+200;
  if (A1.coord_seq == A2.coord_seq) // on s'assure que les 2 reads sont sur la meme sequence
    {
      pos1 = A1.coord_pos;
      pos2 = A2.coord_pos;
      if (pos1 < pos2) // verification des orientations des reads paires
	{
	  if (((pos2 - pos1) > size_insert_min) && ((pos2 - pos1) < size_insert_max)) // verification de l'ecartement entre reads
	    {
	      ret = 1;
	    }
	}
    }
  return ret;
}


// codage des alignements
// le code est recupere dans la chaine de caractere codali
// c'est une suite d'alignements elementaires 
// codage d'une substitution : S pos x   [ pos = position dans le read de la substitution, x = A,C,G ou T ]
// codage d'une insertion :    I pos x+  [ pos = position dans le read de l'insertion, x+ =  1 ou plusieurs caracteres A,C,G ou T ]
// codage d'une deletion :     D pos x+  [ pos = position dans le read de la deletion, x+ =  1 ou plusieurs caracteres A,C,G ou T ]
// le dernier caratere est le symbole X
// example S 12 A D 56 A T G I 87 T C X ==> substitution (A) position 12, deletion (ATG) position 56, insertion (TC) position 87

int code_alignment(int8_t *code, int score, int8_t *gen, int8_t *read)
{
  int i,k,s,x,ix,jx;
  BACKTRACK BT[SIZE_READ];

  // codage substitution          CODE_SUB pos x     
  //        deletion              CODE_DEL pos x+    
  //        insertion             CODE_INS pos x+    
  //        fin                   CODE_END           
  //
  // x  = A | C | G | T    
  // x+ = une sequence d'au moins 1 X (i.e. A, C, G ou T)
  // pos = entier (8 bits) : donne la position du variant par rapport au debut du read
  
  if (score == 0) { code[0]=CODE_END; return 1; }

  // on regarde d'abord si on a a faire a des erreurs de substitution uniquement (c'est vrai la plupart du temps)
  // dans ce cas, on verifie simplement la diagonale de la matrice et on repere les positions de subsitutions
  k = 0;
  s = 0;
  for (i=SIZE_SEED; i<SIZE_NBR4+SIZE_SEED; i++) 
    {
      if (gen[i]!=read[i]) 
	{ 
	  s=s+COST_SUB; 
	  code[k++] = CODE_SUB; 
	  code[k++] = i; 
	  code[k++] = read[i]; 
	  if (s>score) break;  
	} 
    }
  code[k++] = CODE_END;

  if (s==score) return k;


  // sinon, on recalcule la matrice (sur quelques diagonales uniquement)
  // dans BT on recupere le chemin du backtracking
  // x = taille du chemin

  x = DPD(&gen[SIZE_SEED],&read[SIZE_SEED],BT);

  x--; // indice sur le chemin
  k=0; // indice sur le code de l'alignement
  while (x>0)
    {
      if (BT[x].type == CODE_SUB) // substitution
	{
	  code[k++] = CODE_SUB;
	  code[k++] = BT[x].jx+SIZE_SEED-1;
	  code[k++] = read[BT[x].jx+SIZE_SEED-1]; 
	  x--;
	}
      else
	{
	  if (BT[x].type == CODE_DEL) // deletion
	    {
	      code[k++] = CODE_DEL; 
	      code[k++] = BT[x].ix+SIZE_SEED;
	      code[k++] = gen[BT[x].ix+SIZE_SEED];
	      jx = BT[x].jx;
	      x--; 
	      while ((BT[x].type == CODE_DEL) && (jx==BT[x].jx)) 
		{ 
		  code[k++] = gen[BT[x].ix+SIZE_SEED];
		  x--; 
		}
	    }
	  else
	    {
	      code[k++] = CODE_INS; // insertion
	      code[k++] = BT[x].jx+SIZE_SEED-1;
	      code[k++] = read[BT[x].jx+SIZE_SEED];
	      ix = BT[x].ix;
	      x--;
	      while ((BT[x].type == CODE_INS) && (ix==BT[x].ix)) 
		{ 
		  code[k++] = read[BT[x].jx+SIZE_SEED];
		  x--; 
		}
	    }
	}
    }
  code[k++] = CODE_END;
  return k;
}


void set_variant(ALIGN A, GENOME *RG, int8_t *BUF_READ, VARINDEL **INDEL, int *SUB, int8_t *MAPCOV)
{
  int i, ix, type, size, content, snp, ptsnp, x;
  int8_t codali[256];
  VARINDEL *indel, *pindel;

  int pos = RG->pt_seq[A.coord_seq]+A.coord_pos;  // position dans le genome qui correspond au premier charactere du read mappe

  // mise a jour de MAPCOV = nombre de reads qui mapped a cette position dans le genome
  MAPCOV[pos] += 1;

  // codage de l'alignement
  // on recupere dans le tableau codali le code des alignements (voir la fonction code_alignment)
  code_alignment(codali,A.score,&RG->data[pos],&BUF_READ[A.num_read*SIZE_READ]);   

  // mise a jour du champ cntSNP (compteur de SNPs)
  // et gestion des indels
  i = 0;
  while (codali[i]!=CODE_END)
    {
      type = codali[i];
      ptsnp = pos+codali[i+1]; // ptsnp pointe sur la position de l'indels
      if (codali[i]==CODE_SUB)     // variant de type substitution : S pos x  { codali[i], codali[i+1], codali[i+2] }
	{
	  snp = codali[i+2];       // SNP = 0,1,2,3  (code A,C,T,G)
	  x = 1<<(snp*8);          // positionne +1 sur le bon compteur
	  SUB[ptsnp] += x;         // SUB[ptsnp] = entier 32 bits qui contient 4 compteurs 8 bits
	  i = i+3;
	}
      else // variant de type indel : D/I pos x+
        {
          i = i+2;
	  content = 0;
	  size = 0;
          while (codali[i]<4) // lecture des codes des caracteres A, C, G, T
            {
	      content = (content<<2)+codali[i];
              i++; size++;
            }
	  ix = ptsnp & (SIZE_LIST_INDEL -1); 
	  indel = INDEL[ix];
	  while (indel != NULL)
	    {
	      if ((indel->type == type) && (indel->pos_seq == ptsnp) && (indel->num_seq == A.coord_seq) && (indel->content==content))
		{
		  indel->count += 1;
		  break;
		}
	      indel = (VARINDEL *) indel->next;
	    }
	  if (indel == NULL)
	    {
	      pindel = INDEL[ix];
	      INDEL[ix] = (VARINDEL *) malloc(sizeof(VARINDEL));
	      INDEL[ix]->type = type;
	      INDEL[ix]->pos_seq = ptsnp;
	      INDEL[ix]->num_seq = A.coord_seq;
	      INDEL[ix]->size = size;
	      INDEL[ix]->count = 1;
	      INDEL[ix]->content = content;
	      INDEL[ix]->next = pindel;
	    }
        }
    }
}

int process_read(GENOME *RG, int8_t *BUF_READ, int nb_read,  VARINDEL **INDEL, int *SUB, int8_t *MAPCOV, TIMES *CT)
{
    double t1, t2;
    ALIGN *A = (ALIGN *) malloc(sizeof(ALIGN)*MAX_ALIGN*NB_DPU);
    int N,i, k, type, numpair, nb_match, numdpu, num_read, pa, pb;
    long coord;
    int nb_read_map = 0;
    int offset[4], score[4], nbread[4];
    int *TMP = (int *) malloc(sizeof(int)*MAX_ALIGN);

    t1 = my_clock();

    // on recupere dans A les alignements en provenance des DPUs
    nb_match = 0;
    for (numdpu=0; numdpu<NB_DPU; numdpu++)
      {
	k=0;
	while ((num_read = read_out_num(numdpu,k)) != -1)
	  {
	    A[nb_match].num_read = num_read;
	    coord = read_out_coord(numdpu,k);
	    A[nb_match].coord_seq = (int) (coord >> 32);
	    A[nb_match].coord_pos = (int) (coord & 0xFFFFFFFF);
	    A[nb_match].score = read_out_score(numdpu,k);
	    k++;
	    nb_match++;
	  }
      }

    // tri des alignements par leur numero et leur score
    qsort(A,nb_match,sizeof(ALIGN),cmpalign);

    // le numero d'une paire est donnee par num_read / 4 (voir la fonction dispatch)
    // les reads sont types de par leur numero (voir la fonction dispatch)
    // type =  num_read%4 ==  0 ==> read PE1
    //                        1 ==> complement read PE1
    //                        2 == >read PE2
    //                        3 ==> complement read PE2
    // les paires a considerer sont donc les couples [0,3] et [1,2]

    i=0;
    while (i<nb_match)
      {
	// offset[x] = indice dans A du 1er read correspondant au type x
	// nbread[x] = nombre de reads du type x
	// score[x] = score min du type x
	for (k=0; k<4; k++) { score[k] = 1000; offset[k] = -1; nbread[k] = 0; }
	numpair = A[i].num_read / 4;
	while ((i<nb_match) && (numpair == A[i].num_read/4)) // on recupres tous les reads d'une paire
	  {
	    type = A[i].num_read % 4;
	    // normalement, la fonction cmpalign tri par numero de read et par score
	    // ainsi on ne retient que les reads avec les scores min
	    if (A[i].score < score[type]) { score[type] = A[i].score; offset[type] = i; nbread[type] = 1; }
	    else { if (A[i].score == score[type]) nbread[type]++; };
	    i++;
	  }
	N=0;
	// traitement d'une paire [0,3]
	// on regarde la compatibilité de tous les reads de type 0 avec les reads de type 3
	for (pa=offset[0]; pa<offset[0]+nbread[0]; pa++)
	  for(pb=offset[3]; pb<offset[3]+nbread[3]; pb++)
	    if (check_pair1(A[pa],A[pb]))
	      if (N<MAX_ALIGN-2) 
		{ 
		  TMP[N++] = pa; TMP[N++] = pb; 
		}
	// traitement d'une paire [1,2]
	// on regarde la compatibilité de tous les reads de type 1 avec les reads de type 2
	for (pa=offset[1]; pa<offset[1]+nbread[1]; pa++)
	  for(pb=offset[2]; pb<offset[2]+nbread[2]; pb++)
	    if (check_pair2(A[pa],A[pb]))
	      if (N<MAX_ALIGN-2) 
		{ 
		  TMP[N++] = pa; TMP[N++] = pb; 
		}
	// choix drastique (mais qui donne d'excellents résultats!)
	if (N==2)
	  {
	    set_variant(A[TMP[0]],RG,BUF_READ,INDEL,SUB,MAPCOV);
	    set_variant(A[TMP[1]],RG,BUF_READ,INDEL,SUB,MAPCOV);
	    nb_read_map += 2;
	  }
      }

    free(TMP);
    free(A);
    t2 = my_clock();
    CT->process_read = t2-t1;
    CT->tot_process_read += t2-t1;
    return nb_read_map;
}


