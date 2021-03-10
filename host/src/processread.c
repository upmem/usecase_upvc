/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "accumulateread.h"
#include "genome.h"
#include "getread.h"
#include "processread.h"
#include "upvc.h"
#include "vartree.h"

#define SIZE_INSERT_MEAN (400)
#define SIZE_INSERT_STD (3 * 50)

#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13
#define CODE_ERR 14

#define CODE_A 0 /* ('A'>>1)&3   41H  0100 0001 */
#define CODE_C 1 /* ('C'>>1)&3   43H  0100 0011 */
#define CODE_T 2 /* ('T'>>1)&3   54H  0101 0100 */
#define CODE_G 3 /* ('G'>>1)&3   47H  0100 0111 */

#define PQD_INIT_VAL (99)

#define PATH_SUBSTITUTION (0)
#define PATH_INSERTION (1)
#define PATH_DELETION (2)

typedef struct {
    int type;
    int ix;
    int jx;
} backtrack_t;

static int min(int a, int b) { return a < b ? a : b; }

static void DPD_compute(
    int s1, int s2, int *Dij, int Dijm, int Dimj, int Dimjm, int *Pij, int Pijm, int *Qij, int Qimj, int *path)
{
    int min_QP, d;

    *Pij = min(Dijm + COST_GAPO, Pijm + COST_GAPE);
    *Qij = min(Dimj + COST_GAPO, Qimj + COST_GAPE);

    if (*Pij < *Qij) {
        min_QP = *Pij;
        *path = PATH_INSERTION;
    } else {
        min_QP = *Qij;
        *path = PATH_DELETION;
    }
    d = Dimjm;
    if ((s1 & 3) != (s2 & 3)) {
        d += COST_SUB;
    }
    if (d < min_QP) {
        *Dij = d;
        *path = PATH_SUBSTITUTION;
    } else {
        *Dij = min_QP;
    }
}

int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, int size_neighbour_in_symbols)
{
    int matrix_size = size_neighbour_in_symbols + 1;
    int diagonal = (NB_DIAG / 2) + 1;
    int D[matrix_size][matrix_size];
    int P[matrix_size][matrix_size];
    int Q[matrix_size][matrix_size];
    int path[matrix_size][matrix_size];
    int min_score = PQD_INIT_VAL;
    int min_score_i_idx = 0;
    int min_score_j_idx = 0;
    int align_distance = 1;

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            D[i][j] = 0;
        }
    }
    for (int i = 0; i <= diagonal; i++) {
        P[i][0] = PQD_INIT_VAL;
        P[0][i] = PQD_INIT_VAL;
        Q[i][0] = PQD_INIT_VAL;
        Q[0][i] = PQD_INIT_VAL;
        D[i][0] = i * COST_SUB;
        D[0][i] = i * COST_SUB;
    }

    for (int i = 1; i < diagonal; i++) {
        for (int j = 1; j < i + diagonal; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &path[i][j]);
        }
        Q[i][i + diagonal] = PQD_INIT_VAL;
        D[i][i + diagonal] = PQD_INIT_VAL;
    }
    for (int i = diagonal; i < matrix_size - diagonal; i++) {
        P[i][i - diagonal] = PQD_INIT_VAL;
        D[i][i - diagonal] = PQD_INIT_VAL;
        for (int j = i - diagonal + 1; j < i + diagonal; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &path[i][j]);
        }
        Q[i][i + diagonal] = PQD_INIT_VAL;
        D[i][i + diagonal] = PQD_INIT_VAL;
    }

    for (int i = matrix_size - diagonal; i < matrix_size; i++) {
        P[i][i - diagonal] = PQD_INIT_VAL;
        D[i][i - diagonal] = PQD_INIT_VAL;
        for (int j = i - diagonal + 1; j < matrix_size; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &path[i][j]);
        }
        if (D[i][matrix_size - 1] < min_score) {
            min_score = D[i][matrix_size - 1];
            min_score_i_idx = i;
            min_score_j_idx = matrix_size - 1;
        }
    }
    for (int j = matrix_size - diagonal; j < matrix_size; j++) {
        if (D[matrix_size - 1][j] < min_score) {
            min_score = D[matrix_size - 1][j];
            min_score_i_idx = matrix_size - 1;
            min_score_j_idx = j;
        }
    }

    {
        int i = min_score_i_idx;
        int j = min_score_j_idx;

        /* Delete the INDELS at the ends */
        while (path[i][j] != PATH_SUBSTITUTION) {
            if (path[i][j] == PATH_INSERTION) {
                j--;
            } else if (path[i][j] == PATH_DELETION) {
                i--;
            } else {
                return -1;
            }
        }

        /* i>1 && j>1 conditions erased the INDELS at the beginning */
        backtrack[0].type = CODE_END;
        while ((i > 1) && (j > 1)) {
            if (path[i][j] == PATH_SUBSTITUTION) {
                i--;
                j--;
                if (D[i][j] != D[i - 1][j - 1]) {
                    backtrack[align_distance].type = CODE_SUB;
                    backtrack[align_distance].ix = i;
                    backtrack[align_distance].jx = j;
                    align_distance++;
                }
            } else {
                if (path[i][j] == PATH_INSERTION) {
                    j--;
                    backtrack[align_distance].type = CODE_INS;
                    backtrack[align_distance].ix = i;
                    backtrack[align_distance].jx = j;
                    align_distance++;
                } else if (path[i][j] == PATH_DELETION) {
                    i--;
                    backtrack[align_distance].type = CODE_DEL;
                    backtrack[align_distance].ix = i;
                    backtrack[align_distance].jx = j;
                    align_distance++;
                } else {
                    ERROR_EXIT(ERR_PROCESSREAD_DPD_FAILED, "Error during DPD compute");
                }
            }
        }
    }

    return align_distance;
}

/*
 * encoding of differences between read and sequence of the reference genome.
 * coding: substitution          CODE_SUB pos x
 *         deletion              CODE_DEL pos x+
 *         insertion             CODE_INS pos x+
 *         end                   CODE_END
 *
 * x  = A | C | G | T
 * x+ = a sequence of at least 1 element (i.e. A, C, G ou T)
 * pos = integer (8 bits) : give the offset of the variant from the start of the read
 *
 * example S 12 A D 56 A T G I 87 T C X ==> substitution (A) position 12, deletion (ATG) position 56, insertion (TC) position 87
 * The code is return in "code" as a table of int8_t
 */

#if 0
static int code_alignment(uint8_t *code, int score, int8_t *gen, int8_t *read, unsigned size_neighbour_in_symbols)
{
    int code_idx, computed_score, backtrack_idx;
    int size_read = SIZE_READ;
    int size_neighbour = size_neighbour_in_symbols;
    backtrack_t backtrak[size_read];

    if (score == 0) {
        code[0] = CODE_END;
        return 1;
    }

    /* First, looking for subsititution only */
    code_idx = 0;
    computed_score = 0;
    for (int i = SIZE_SEED; i < size_neighbour + SIZE_SEED; i++) {
        if ((gen[i] & 3) != read[i]) {
            computed_score += COST_SUB;
            code[code_idx++] = CODE_SUB;
            code[code_idx++] = i;
            code[code_idx++] = read[i];
            if (computed_score > score) {
                break;
            }
        }
    }
    code[code_idx++] = CODE_END;
    if (computed_score == score)
        return code_idx;

    /* Otherwise, re-compute the matrix (only some diagonals) and put in backtrack the path */
    backtrack_idx = DPD(&gen[SIZE_SEED], &read[SIZE_SEED], backtrak, size_neighbour_in_symbols);
    if (backtrack_idx == -1) {
        code[0] = CODE_ERR;
        return 1;
    }

    backtrack_idx--;
    code_idx = 0;
    while (backtrack_idx > 0) {
        if (backtrak[backtrack_idx].type == CODE_SUB) {
            code[code_idx++] = CODE_SUB;
            code[code_idx++] = backtrak[backtrack_idx].jx + SIZE_SEED - 1;
            code[code_idx++] = read[backtrak[backtrack_idx].jx + SIZE_SEED - 1];
            backtrack_idx--;
        } else {
            if (backtrak[backtrack_idx].type == CODE_DEL) {
                int backtrack_jx = backtrak[backtrack_idx].jx;
                code[code_idx++] = CODE_DEL;
                code[code_idx++] = backtrak[backtrack_idx].ix + SIZE_SEED;
                code[code_idx++] = gen[backtrak[backtrack_idx].ix + SIZE_SEED] & 3;
                backtrack_idx--;
                while ((backtrak[backtrack_idx].type == CODE_DEL) && (backtrack_jx == backtrak[backtrack_idx].jx)) {
                    code[code_idx++] = gen[backtrak[backtrack_idx].ix + SIZE_SEED] & 3;
                    backtrack_idx--;
                }
            } else {
                int backtrack_ix = backtrak[backtrack_idx].ix;
                code[code_idx++] = CODE_INS;
                code[code_idx++] = backtrak[backtrack_idx].jx + SIZE_SEED - 1;
                code[code_idx++] = read[backtrak[backtrack_idx].jx + SIZE_SEED];
                backtrack_idx--;
                while ((backtrak[backtrack_idx].type == CODE_INS) && (backtrack_ix == backtrak[backtrack_idx].ix)) {
                    code[code_idx++] = read[backtrak[backtrack_idx].jx + SIZE_SEED];
                    backtrack_idx--;
                }
            }
        }
    }
    code[code_idx++] = CODE_END;
    return code_idx;
}

static void set_variant(
    dpu_result_out_t result_match, genome_t *ref_genome, int8_t *reads_buffer, unsigned int size_neighbour_in_symbols)
{
    uint32_t code_result_idx;
    uint8_t code_result_tab[256];
    int8_t *read;
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    uint64_t genome_pos = ref_genome->pt_seq[result_match.coord.seq_nr] + result_match.coord.seed_nr;
    int size_read = SIZE_READ;

    /* Get the differences betweend the read and the sequence of the reference genome that match */
    read = &reads_buffer[result_match.num * size_read];
    code_alignment(code_result_tab, result_match.score, &ref_genome->data[genome_pos], read, size_neighbour_in_symbols);
    if (code_result_tab[0] == CODE_ERR)
        return;

    /* Update "mapping_coverage" with the number of reads that match at this position of the genome */
    for (int i = 0; i < size_read; i++) {
        ref_genome->mapping_coverage[genome_pos + i] += 1;
    }

    code_result_idx = 0;
    while (code_result_tab[code_result_idx] != CODE_END) {
        int code_result = code_result_tab[code_result_idx];
        int64_t pos_variant_read = code_result_tab[code_result_idx + 1];
        int64_t pos_variant_genome = genome_pos + pos_variant_read;
        int ref_pos = 0;
        int alt_pos = 0;
        variant_t *newvar = (variant_t *)malloc(sizeof(variant_t));
        newvar->depth = 1;
        newvar->score = result_match.score;
        newvar->next = NULL;
        if (code_result == CODE_SUB) {
            /* SNP = 0,1,2,3  (code A,C,T,G) */
            int snp = code_result_tab[code_result_idx + 2];
            newvar->ref[ref_pos++] = nucleotide[ref_genome->data[pos_variant_genome] & 3];
            newvar->alt[alt_pos++] = nucleotide[snp & 3];

            code_result_idx += 3;
        } else if (code_result == CODE_INS) {
            int64_t ps_var_genome = pos_variant_genome;
            int64_t ps_var_read = pos_variant_read;
            code_result_idx += 2;

            while (code_result_tab[code_result_idx] < 4) {
                ps_var_read++;
                code_result_idx++;
            }

            while (ref_genome->data[ps_var_genome] == read[ps_var_read]) {
                ps_var_genome--;
                ps_var_read--;
                pos_variant_genome--;
                pos_variant_read--;
            }

            newvar->ref[ref_pos++] = nucleotide[ref_genome->data[pos_variant_genome] & 3];

            while (pos_variant_read <= ps_var_read) {
                newvar->alt[alt_pos++] = nucleotide[read[pos_variant_read] & 3];
                if (alt_pos >= MAX_SIZE_ALLELE - 1) {
                    free(newvar);
                    return;
                }
                pos_variant_read++;
            }

        } else if (code_result == CODE_DEL) {
            int64_t ps_var_genome = pos_variant_genome;
            int64_t ps_var_read = pos_variant_read;
            code_result_idx += 2;

            while (code_result_tab[code_result_idx] < 4) {
                ps_var_genome++;
                code_result_idx++;
            }

            while (ref_genome->data[ps_var_genome] == read[ps_var_read]) {
                ps_var_read--;
                ps_var_genome--;
                pos_variant_genome--;
                pos_variant_read--;
            }

            newvar->alt[alt_pos++] = nucleotide[ref_genome->data[pos_variant_genome] & 3];

            while (pos_variant_genome <= ps_var_genome) {
                newvar->ref[ref_pos++] = nucleotide[ref_genome->data[pos_variant_genome] & 3];
                if (ref_pos >= MAX_SIZE_ALLELE - 1) {
                    free(newvar);
                    return;
                }
                pos_variant_genome++;
            }
            pos_variant_genome -= ref_pos;
        }
        newvar->ref[ref_pos] = '\0';
        newvar->alt[alt_pos] = '\0';
        variant_tree_insert(
            newvar, result_match.coord.seq_nr, pos_variant_genome + 1 - ref_genome->pt_seq[result_match.coord.seq_nr]);
    }
}
#endif

static pthread_mutex_t non_mapped_mutex;
static void add_to_non_mapped_read(int numread, int round, FILE *fpe1, FILE *fpe2, int8_t *reads_buffer)
{
    if (fpe1 == NULL || fpe2 == NULL)
        return;
    pthread_mutex_lock(&non_mapped_mutex);
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    int size_read = SIZE_READ;
    int8_t *read = &reads_buffer[numread * size_read];
    fprintf(fpe1, ">>%d\n", SIZE_SEED * (round + 1));
    for (int j = SIZE_SEED; j < size_read; j++) {
        fprintf(fpe1, "%c", nucleotide[read[j] & 3]);
    }
    for (int j = 0; j < SIZE_SEED; j++) {
        fprintf(fpe1, "A");
    }
    fprintf(fpe1, "\n");
    read = &reads_buffer[(numread + 2) * size_read];
    fprintf(fpe2, ">>%d\n", SIZE_SEED * (round + 1));
    for (int j = SIZE_SEED; j < size_read; j++) {
        fprintf(fpe2, "%c", nucleotide[read[j] & 3]);
    }
    for (int j = 0; j < SIZE_SEED; j++) {
        fprintf(fpe2, "A");
    }
    fprintf(fpe2, "\n");
    pthread_mutex_unlock(&non_mapped_mutex);
}

static void update_frequency_table(
    genome_t *ref_genome, 
    dpu_result_out_t *result_tab, 
    int8_t *reads_buffer, 
    int pos) {

  struct frequency_info **frequency_table = get_frequency_table();
  uint64_t genome_pos = ref_genome->pt_seq[result_tab[pos].coord.seq_nr] + result_tab[pos].coord.seed_nr; 
  int num = result_tab[pos].num;
  int8_t *read = reads_buffer + (num * SIZE_READ);
  for(int j = 0; j < SIZE_READ; ++j) {
    if(genome_pos + j < genome_get()->fasta_file_size) {
      frequency_table[read[j]][genome_pos+j].freq++;
      frequency_table[read[j]][genome_pos+j].score += result_tab[pos].score;
    }
    else
      printf("WARNING: reads matched at position that exceeds genome size\n");
  }
}

static volatile unsigned int curr_match;
static pthread_mutex_t curr_match_mutex;
unsigned int acquire_curr_match()
{
    pthread_mutex_lock(&curr_match_mutex);
    return curr_match;
}
void release_curr_match(unsigned int new_curr_match)
{
    curr_match = new_curr_match;
    pthread_mutex_unlock(&curr_match_mutex);
    return;
}

static pthread_barrier_t barrier;

typedef struct {
    unsigned int nb_match;
    dpu_result_out_t *result_tab;
    int round;
    int8_t *reads_buffer;
    genome_t *ref_genome;
    FILE *fpe1;
    FILE *fpe2;
} process_read_arg_t;

static uint64_t nr_reads_total = 0ULL;
static uint64_t nr_reads_non_mapped = 0ULL;
static pthread_mutex_t nr_reads_mutex = PTHREAD_MUTEX_INITIALIZER;

static void do_process_read(process_read_arg_t *arg)
{
    const unsigned int nb_match = arg->nb_match;
    dpu_result_out_t *result_tab = arg->result_tab;
    int round = arg->round;
    int8_t *reads_buffer = arg->reads_buffer;
    genome_t *ref_genome = arg->ref_genome;
    FILE *fpe1 = arg->fpe1;
    FILE *fpe2 = arg->fpe2;
    //unsigned int size_neighbour_in_symbols = (SIZE_NEIGHBOUR_IN_BYTES - DELTA_NEIGHBOUR(round)) * 4;

    /*
     * The number of a pair is given by "num_read / 4 " (see dispatch_read function)
     * Their type is given by their offset (see dispatch_read function)
     * type = num_read%4 == 0 ==> read 1
     *                      1 ==> read 1 complement
     *                      2 ==> read 2
     *                      3 ==> read 2 complement
     * The read pair to consider are [0, 3] and [1, 2].
     *
     * NEW (04/2020):
     *  - more paired reads are considered
     *  - when different position mapping are possible, choose the less covered zone
     */

    while (true) {
        unsigned int i;
        if ((i = acquire_curr_match()) >= nb_match) {
            release_curr_match(i);
            return;
        }
        int numpair = result_tab[i].num / 4;
        unsigned int j = i;
        while ((j < nb_match) && (numpair == result_tab[j].num / 4)) {
            j++;
        }
        release_curr_match(j);

        // i = start index in result_tab
        // j = stop index in result_tab
        // select best couples of paired reads
        unsigned int P1[1000];
        unsigned int P2[1000];
        unsigned int np = 0;
        unsigned int P1_all[1000];
        unsigned int P2_all[1000];
        unsigned int np_all = 0;
        unsigned int pos1, pos2, t1, t2;
        unsigned int best_score = 1000;
        unsigned int best_score_all = 1000;
        // test all significant pairs of reads (0,3) & (1,2)
        for (unsigned int x1 = i; x1 < j; x1++) {
            t1 = result_tab[x1].num % 4;
            pos1 = result_tab[x1].coord.seed_nr;
            for (unsigned int x2 = i + 1; x2 < j; x2++) {
                pos2 = result_tab[x2].coord.seed_nr;
                t2 = result_tab[x2].num % 4;
                if (t1 + t2 == 3) // select significant pair
                {
                    if ((abs((int)pos2 - (int)pos1) > 50 && (abs((int)pos2 - (int)pos1) < 1000))) {
                        if (result_tab[x1].score + result_tab[x2].score < best_score) {
                            np = 0;
                            best_score = result_tab[x1].score + result_tab[x2].score;
                            P1[np] = x1;
                            P2[np] = x2;
                            np++;
                        } else {
                            if (result_tab[x1].score + result_tab[x2].score == best_score) {
                                P1[np] = x1;
                                P2[np] = x2;
                                if (np < 999)
                                    np++;
                            }
                        }
                    }
                    else {
                      
                        if (result_tab[x1].score + result_tab[x2].score < best_score_all) {
                            np_all = 0;
                            best_score_all = result_tab[x1].score + result_tab[x2].score;
                            P1_all[np_all] = x1;
                            P2_all[np_all] = x2;
                            np_all++;
                        } else if (result_tab[x1].score + result_tab[x2].score == best_score_all) {
                                P1_all[np_all] = x1;
                                P2_all[np_all] = x2;
                                if (np_all < 999)
                                    np_all++;
                        }
                    }
                }
            }
        }
        if (np > 0) {

            // update frequency table 
            for (unsigned int i = 0; i < np; i++) {

              // for each matching position of this read in the reference genome, add +1 in the corresponding nucleotide column
              update_frequency_table(ref_genome, result_tab, reads_buffer, P1[i]);
              update_frequency_table(ref_genome, result_tab, reads_buffer, P2[i]);
            }
        }
        else if (np_all > 0) {

            // update frequency table 
            for (unsigned int i = 0; i < np_all; i++) {

              // for each matching position of this read in the reference genome, add +1 in the corresponding nucleotide column
              update_frequency_table(ref_genome, result_tab, reads_buffer, P1_all[i]);
              update_frequency_table(ref_genome, result_tab, reads_buffer, P2_all[i]);
            }
        } else {
            pthread_mutex_lock(&nr_reads_mutex);
            nr_reads_non_mapped++;
            pthread_mutex_unlock(&nr_reads_mutex);
            add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
        }
        pthread_mutex_lock(&nr_reads_mutex);
        nr_reads_total++;
        pthread_mutex_unlock(&nr_reads_mutex);
    }
}

#define PROCESS_READ_THREAD (8)
#define PROCESS_READ_THREAD_SLAVE (PROCESS_READ_THREAD - 1)
static process_read_arg_t args;
static pthread_t thread_id[PROCESS_READ_THREAD_SLAVE];
static bool stop_threads = false;

void process_read(FILE *fpe1, FILE *fpe2, int round, unsigned int pass_id)
{
    int8_t *reads_buffer = get_reads_buffer(pass_id);
    acc_results_t acc_res = accumulate_get_result(pass_id);

    curr_match = 0;

    args.nb_match = acc_res.nb_res;
    args.result_tab = acc_res.results;
    args.round = round;
    args.reads_buffer = reads_buffer;
    args.fpe1 = fpe1;
    args.fpe2 = fpe2;

    pthread_barrier_wait(&barrier);
    do_process_read(&args);
    pthread_barrier_wait(&barrier);

    free(acc_res.results);
}

static void *process_read_thread_fct(void *arg)
{
    pthread_barrier_wait(&barrier);
    while (!stop_threads) {
        do_process_read(arg);
        pthread_barrier_wait(&barrier);
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

void process_read_init()
{
    genome_t *ref_genome = genome_get();
    args.ref_genome = ref_genome;

    assert(pthread_mutex_init(&curr_match_mutex, NULL) == 0);
    assert(pthread_mutex_init(&non_mapped_mutex, NULL) == 0);
    assert(pthread_barrier_init(&barrier, NULL, PROCESS_READ_THREAD) == 0);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_create(&thread_id[each_thread], NULL, process_read_thread_fct, &args) == 0);
    }
}

void process_read_free()
{
    stop_threads = true;
    pthread_barrier_wait(&barrier);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_join(thread_id[each_thread], NULL) == 0);
    }

    assert(pthread_barrier_destroy(&barrier) == 0);
    assert(pthread_mutex_destroy(&curr_match_mutex) == 0);
    assert(pthread_mutex_destroy(&non_mapped_mutex) == 0);
    fprintf(stderr, "%% reads non mapped: %f%%\n", (float)nr_reads_non_mapped * 100.0 / (float)nr_reads_total);
}
