#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

#include "upvc.h"

int SIZE_READ;
int SIZE_NBR;
int SIZE_NBR4;

inline double my_clock(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (1.0e-6*t.tv_usec + t.tv_sec);
}

void static *align_on_dpu (void *a)
{
  ARG *arg = (ARG *) a;
  align(arg->arg1);
  return NULL;
}

void run_dpu(char *prog, TIMES *CT)
{
  int numdpu;
  double t1, t2;
  pthread_t THREAD_DPU[NB_DPU];                 // threads pour paralleliser le calcul sur les DPUs
  ARG *arg_align[NB_DPU];                       // structure pour passer les arguments                                       

  t1 = my_clock();

  for (numdpu=0; numdpu<NB_DPU; numdpu++)
    {
      arg_align[numdpu] = (ARG *) malloc(sizeof(ARG));
      arg_align[numdpu]->arg1 = numdpu;
    }

  if (strcmp(prog,"mapping") == 0)
    {
      for (numdpu=0; numdpu<NB_DPU; numdpu++)
        {
          pthread_create(&THREAD_DPU[numdpu],NULL,align_on_dpu,arg_align[numdpu]);
        }
      for (numdpu=0; numdpu<NB_DPU; numdpu++)
        {
          pthread_join(THREAD_DPU[numdpu],NULL);
        }
    }
  for (numdpu=0; numdpu<NB_DPU; numdpu++) free(arg_align[numdpu]);
  t2 = my_clock();
  CT->map_read = t2-t1;
  CT->tot_map_read += t2-t1;
}

int main (int argc, char *argv[])
{
  char filename[1024];
  GENOME *REFGENOME;               // genome de reference
  int8_t *MAPCOV;                  // couverture du mapping
  INDEX_SEED **SEED;               // index des graines
  VARINDEL **INDEL;                // liste des insertions/deletions
  int *SUB;                        // liste des substitutions
  TIMES *CT;                       // temps de calcul

  if (argc != 2) { printf ("\nusage:\n  %s <genome> \n\n",argv[0]); exit (255); }

  // initialisation des compteurs de temps de calcul
  CT = (TIMES *) malloc(sizeof(TIMES));
  CT->get_genome = 0.0;
  CT->index_genome = 0.0;
  CT->get_reads = 0.0;
  CT->tot_get_reads = 0.0;
  CT->dispatch_read = 0.0;
  CT->tot_dispatch_read = 0.0;
  CT->map_read = 0.0;
  CT->tot_map_read = 0.0;
  CT->process_read = 0.0;
  CT->tot_process_read = 0.0;
  CT->vcf = 0.0;  

  // initialisation des variables externes
  SIZE_READ = getSizeRead(argv[1]);        // lecture 1er read d'un fichier PE. Hypothese : tous les reads sont de meme taille - getread.c
  SIZE_NBR = (SIZE_READ-SIZE_SEED)/4;      // taille du voisinage en octet
  SIZE_NBR4 = SIZE_NBR*4;

  printf ("Information\n");
  printf (" - read size: %d\n",SIZE_READ);

  //  lecture du genome
  printf ("Read genome\n");
  REFGENOME = get_genome(argv[1],CT);  // getgenome.c
  printf (" - #seq: %d\n",REFGENOME->nb_seq);
  printf (" - time: %lf\n",CT->get_genome);

  // indexation du genome
  printf ("Index genome\n");
  SEED = index_genome(REFGENOME,CT); // index.c
  printf (" - time: %lf\n",CT->index_genome);

  // initialisation de la couverture du mapping
  MAPCOV = (int8_t *) calloc(sizeof(int8_t), REFGENOME->sizefile);

  // initialisation (a zero) des compteurs de substitution
  SUB = (int *) calloc(sizeof(int), REFGENOME->sizefile);

  // initialisation liste des insertions/omissions
  INDEL = (VARINDEL **) malloc (sizeof(VARINDEL *)*SIZE_LIST_INDEL);
  { 
    int i;
    for (i=0; i<SIZE_LIST_INDEL; i++) INDEL[i] = NULL;
  }
  // boucle principale
  //   - lecture d'un packet de reads
  //   - dispatch des reads de ce packet dans les DPUs
  //   - execution du mapping sur les DPUs
  //   - postprocessing des alignements

  {
    FILE *fpe1, *fpe2;                          // descripteurs de fichier pair-end 1 & 2
    int8_t *BUF_READ;                           // buffer pour stocker les reads
    int  nb_read;
    int  nb_read_total = 0;
    int  nb_read_map = 0;
    int  nb_pass = 0;

    // ouverture des fichiers de reads
    sprintf(filename,"%s_PE1.fastq",argv[1]);
    fpe1 = fopen(filename,"r");
    sprintf(filename,"%s_PE2.fastq",argv[1]);
    fpe2 = fopen(filename,"r");

    // initialisation du buffer de reads
    BUF_READ = (int8_t *) malloc(sizeof(int8_t)*MAX_BUF_READ*SIZE_READ);

    // lecture packet de reads
    while ((nb_read = getPEreads(fpe1, fpe2, BUF_READ,CT)) != 0) // getread.c
      {
	printf ("Pass %d\n",nb_pass);

	nb_read_total += nb_read;
	printf (" - get %d reads (%d)\n",nb_read/2,nb_read_total/2);
	printf (" - time to get reads      : %7.2lf sec. / %7.2lf sec.\n",CT->get_reads,CT->tot_get_reads);

	// dispatch des reads et de leur complement dans les DPUs
	dispatch_read(SEED, BUF_READ, nb_read, CT);  // dispatch.c
	printf (" - time to dispatch reads : %7.2lf sec. / %7.2lf sec.\n",CT->dispatch_read,CT->tot_dispatch_read);

	// execution du mapping sur les DPUs
	run_dpu("mapping",CT);
	printf (" - time to map reads      : %7.2lf sec. / %7.2lf sec.\n",CT->map_read,CT->tot_map_read);

	// post processing des reads
	nb_read_map += process_read(REFGENOME,BUF_READ,nb_read,INDEL,SUB,MAPCOV,CT); // processread.c
	printf (" - time to process reads  : %7.2lf sec. / %7.2lf sec.\n",CT->process_read,CT->tot_process_read);
	printf (" - map %d reads\n",nb_read_map);

	printf ("\n");
	nb_pass++;
      }
    free(BUF_READ);	
    fclose(fpe1);
    fclose(fpe2);
  }


  // creation fichier vcf
  {
    int nb_var;
    sprintf(filename,"%s.vcf",argv[1]);
    printf ("Create VCF\n");
    nb_var = create_vcf(filename,REFGENOME,INDEL,SUB,MAPCOV,CT);
    printf (" - number of variants: %d\n",nb_var);
    printf (" - time: %lf sec.\n",CT->vcf);
  }


  // restitution de la memoire
  free_genome(REFGENOME);
  free_index(SEED);
  free(SUB);
  free(CT);
  free_dpu();
  return 0;
}
