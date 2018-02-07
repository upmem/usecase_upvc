#define VERSION           "VERSION 1.3"
#define SIZE_SEED         12
#define NB_SEED           16777216
#define NB_DPU            128
#define SIZE_IDX_HOST     (NB_SEED/NB_DPU)
#define MAX_NB_DPU_READ   16384                  // nombre de reads par DPU par passe
#define MAX_BUF_READ      1048576                // nombre de reads total par passe
#define MAX_ALIGN         65536
#define MAX_SEQ_GEN       1000
#define MAX_SIZE_IDX_SEED 1000
#define SIZE_WINDOW_CPLX  11

#define SIZE_LIST_INDEL  (1<<16)

#define COST_SUB         10
#define COST_GAPO        11
#define COST_GAPE        1
#define MAX_SCORE        40
#define NB_DIAG          15

#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13
#define CODE_A    0    // ('A'>>1)&3   41H  0100 0001
#define CODE_C    1    // ('C'>>1)&3   43H  0100 0011
#define CODE_T    2    // ('T'>>1)&3   54H  0101 0100
#define CODE_G    3    // ('G'>>1)&3   47H  0100 0111

extern int SIZE_NBR;
extern int SIZE_NBR4;
extern int SIZE_READ;

typedef struct {
  int arg1;
  int arg2;
} ARG;

// structure de donnees associee au genome de reference
// initialisee dans get_genome() - getgenome.c
typedef struct {
  int8_t    *data;                 // Chaine d'entiers (8 bits) qui contient le genome de reference
  int        nb_seq;               // Nombre de sequences qui composent le genome
  long      *pt_seq;               // pointeur dans data sur le debut des sequences
  int       *len_seq;              // taille des séquences
  long      sizefile;              // taille du fichier fasta
} GENOME;

// index associe a une graine
typedef struct INDEX_SEED INDEX_SEED;
struct INDEX_SEED {
  int nb_nbr;        // nombre de voisinage
  int offset;        // adresse dans la memoire DPU du premier voisinage a traiter
  int num_dpu;       // numero du DPU ou est dispatche l'index
  INDEX_SEED *next;  // liste chaine d'index
};


// memorisation des differents temps de calcul
typedef struct {
  double get_genome;
  double index_genome;
  double get_reads;
  double dispatch_read;
  double map_read;
  double process_read;
  double tot_get_reads;
  double tot_dispatch_read;
  double tot_map_read;
  double tot_process_read;
  double vcf;
} TIMES;

// alignement
typedef struct {
  int num_read;
  int coord_seq;
  int coord_pos;
  int score;
} ALIGN;

// backtracking
typedef struct {
  int type;
  int ix;
  int jx;
} BACKTRACK;

// variant de type indel
typedef struct {
  int num_seq;   // numero de sequence
  int pos_seq;   // position du variant sur la sequence
  int type;      // code du type de variant : CODE_SUB, CODE_DEL, CODE_INS
  int size;      // taille du variant
  int content;   // chaine de ATGC codee (max 16)x
  int count;     // nombre d'occurences du variant
  void *next;    // pointeur sur un autre variant
} VARINDEL;

double my_clock(void);

GENOME *get_genome(char* name, TIMES *CT);  // getgenome.c
void free_genome(GENOME *G);  

INDEX_SEED** index_genome(GENOME *G, TIMES *CT); // index.c
void free_index(INDEX_SEED **SEED);

int  create_vcf(char *filename, GENOME *RG, VARINDEL **INDEL, int *SUB, int8_t *MAPCOV, TIMES *CT); // vcf.c

int  getPEreads(FILE *fpe1, FILE *fpe2, int8_t *BUF_READ, TIMES *CT);
int  getSeqFastQ(FILE *ff, int8_t *read1, int8_t *read2);
void dispatch_read(INDEX_SEED **SEED, int8_t* BUF_READ, int nb_pair, TIMES *CT);
int  process_read(GENOME *RG, int8_t *BUF_READ, int nb_read,  VARINDEL **INDEL, int *SUB, int8_t* MAPCOV, TIMES *CT);
int  getSizeRead (char *name);
void align(int d);

int  DPD(int8_t *s1, int8_t *s2, BACKTRACK *bt);

int  code_seed(int8_t *SEQ);
void code_neighbor(int8_t *SEQ, int8_t *CODE);
void decode_neighbor(int8_t *CODE, int8_t *SEQ);

void malloc_neighbor_idx (int d, int k);
void malloc_dpu();
void free_dpu();

// les fonctions suivantes emulent les transactions de/vers les DPU
// elles sont definies dans upvc_dpu.c

// presentes dans index.c
void write_neighbor_idx  (int d, int k, int8_t *v);
void write_coordinate    (int d, int k, long v);   

// presentent dans dispatch.c
void write_neighbor_read (int d, int k, int8_t *v);
void write_count (int d, int k, int v);
void write_offset(int d, int k, int v);
void write_num   (int d, int k, int v);

// presentent dans processread.c
int  read_out_num (int d, int k);
long read_out_coord (int d, int k);
int  read_out_score (int d, int k);


