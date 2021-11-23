/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "accumulateread.h"
#include "common.h"
#include "debug.h"
#include "genome.h"
#include "getread.h"
#include "processread.h"
#include "sam.h"
#include "upvc.h"
#include "vartree.h"

#define DEBUG_READ_MAPPING true

#define SIZE_INSERT_MEAN (400)
#define SIZE_INSERT_STD (3 * 50)

#define CODE_A 0 /* ('A'>>1)&3   41H  0100 0001 */
#define CODE_C 1 /* ('C'>>1)&3   43H  0100 0011 */
#define CODE_T 2 /* ('T'>>1)&3   54H  0101 0100 */
#define CODE_G 3 /* ('G'>>1)&3   47H  0100 0111 */

#define PQD_INIT_VAL (99)

#define MAX_SUBSTITUTION (4)

static bool flag_dbg = false;

typedef struct {
    int type;
    int ix;
    int jx;
} backtrack_t;

static int min(int a, int b) { return a < b ? a : b; }

static void DPD_compute(
    int s1, int s2, int *Dij, int Dijm, int Dimj, int Dimjm, int *Pij, int Pijm, int *Qij, int Qimj, int *xij)
{
    int min_QP, d;

    *Pij = min(Dijm + COST_GAPO, Pijm + COST_GAPE);
    *Qij = min(Dimj + COST_GAPO, Qimj + COST_GAPE);
    *xij = 0;

    int x;
    if (*Pij < *Qij) {
        min_QP = *Pij;
        x = 2;
    } else {
        min_QP = *Qij;
        x = 3;
    }
    d = Dimjm;
    if ((s1 & 3) != (s2 & 3)) {
        d += COST_SUB;
        *xij = 1;
    }
    if (d < min_QP) {
        *Dij = d;
    } else {
        *Dij = min_QP;
        *xij = x;
    }
}

int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, int size_neighbour_in_symbols)
{
    int matrix_size = size_neighbour_in_symbols + 1;
    int diagonal = (NB_DIAG / 2) + 1;
    int D[matrix_size][matrix_size];
    int P[matrix_size][matrix_size];
    int Q[matrix_size][matrix_size];
    int X[matrix_size][matrix_size];				
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
                &Q[i][j], Q[i - 1][j], &X[i][j]);
        }
        Q[i][i + diagonal] = PQD_INIT_VAL;
        D[i][i + diagonal] = PQD_INIT_VAL;
    }
    for (int i = diagonal; i < matrix_size - diagonal; i++) {
        P[i][i - diagonal] = PQD_INIT_VAL;
        D[i][i - diagonal] = PQD_INIT_VAL;
        for (int j = i - diagonal + 1; j < i + diagonal; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &X[i][j]);
        }
        Q[i][i + diagonal] = PQD_INIT_VAL;
        D[i][i + diagonal] = PQD_INIT_VAL;
    }

    for (int i = matrix_size - diagonal; i < matrix_size; i++) {
        P[i][i - diagonal] = PQD_INIT_VAL;
        D[i][i - diagonal] = PQD_INIT_VAL;
        for (int j = i - diagonal + 1; j < matrix_size; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &X[i][j]);
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
        backtrack[0].type = CODE_END;
        while ((i > 0) && (j > 0)) {
            if(X[i][j] == 0) {
                i--;
                j--;
            } else {
                if(X[i][j] == 1) {
                    i--;
                    j--;
                    backtrack[align_distance].type = CODE_SUB;
                    backtrack[align_distance].ix = i;
                    backtrack[align_distance].jx = j;
                    align_distance++;
                } else {
                    if(X[i][j] == 2) {
                        j--;
                        backtrack[align_distance].type = CODE_INS;
                        backtrack[align_distance].ix = i;
                        backtrack[align_distance].jx = j;
                        align_distance++;
                    }
                    else {
                        i--;
                        backtrack[align_distance].type = CODE_DEL;
                        backtrack[align_distance].ix = i;
                        backtrack[align_distance].jx = j;
                        align_distance++;
                    }
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

#ifdef USE_INDEL
static int code_alignment(uint8_t *code, int score, int8_t *gen, int8_t *read, unsigned size_neighbour_in_symbols, bool *flag)
{
    int code_idx, computed_score, backtrack_idx;
    int size_read = SIZE_READ;
    int size_neighbour = size_neighbour_in_symbols;
    backtrack_t backtrak[size_read];

    *flag = false;

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
    backtrack_idx = DPD(gen, read, backtrak, size_neighbour_in_symbols + SIZE_SEED);
    if (backtrack_idx == -1) {
        code[0] = CODE_ERR;
        return 1;
    }

    backtrack_idx--;
    code_idx = 0;
    while (backtrack_idx > 0) {
        if (backtrak[backtrack_idx].type == CODE_SUB) {
            code[code_idx++] = CODE_SUB;
            code[code_idx++] = backtrak[backtrack_idx].jx - 1;
            code[code_idx++] = read[backtrak[backtrack_idx].jx - 1];
            backtrack_idx--;
        } else {
            if (backtrak[backtrack_idx].type == CODE_DEL) {
                //int backtrack_jx = backtrak[backtrack_idx].jx;
                code[code_idx++] = CODE_DEL;
                code[code_idx++] = backtrak[backtrack_idx].ix;
                code[code_idx++] = gen[backtrak[backtrack_idx].ix] & 3;
                backtrack_idx--;
                /*
                while ((backtrak[backtrack_idx].type == CODE_DEL) && (backtrack_jx == backtrak[backtrack_idx].jx)) {
                    code[code_idx++] = gen[backtrak[backtrack_idx].ix] & 3;
                    backtrack_idx--;
                }
                */
            } else {
                //int backtrack_ix = backtrak[backtrack_idx].ix;
                code[code_idx++] = CODE_INS;
                code[code_idx++] = backtrak[backtrack_idx].jx - 1;
                code[code_idx++] = read[backtrak[backtrack_idx].jx];
                backtrack_idx--;
                /*
                while ((backtrak[backtrack_idx].type == CODE_INS) && (backtrack_ix == backtrak[backtrack_idx].ix)) {
                    code[code_idx++] = read[backtrak[backtrack_idx].jx];
                    backtrack_idx--;
                }
                */
            }
        }
    }
    code[code_idx++] = CODE_END;
    return code_idx;
}
#endif

#if 0
static void set_variant(
    dpu_result_out_t result_match, genome_t *ref_genome, int8_t *reads_buffer, unsigned int size_neighbour_in_symbols)
{
    uint32_t code_result_idx;
    uint8_t code_result_tab[256];
    int8_t *read;
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    uint64_t genome_pos = ref_genome->pt_seq[result_match.coord.seq_nr] + result_match.coord.seed_nr;
    int size_read = SIZE_READ;
    //LOG_TRACE("set_variant called\n");

    /* Get the differences betweend the read and the sequence of the reference genome that match */
    read = &reads_buffer[result_match.num * size_read];
    code_alignment(code_result_tab, result_match.score, &ref_genome->data[genome_pos], read, size_neighbour_in_symbols);
    if (code_result_tab[0] == CODE_ERR)
        return;

    /* Update "mapping_coverage" with the number of reads that match at this position of the genome */
    for (int i = 0; i < size_read; i++) {
        ref_genome->mapping_coverage[genome_pos + i] += 1;
    }

#if DEBUG_READ_MAPPING
    // TODO: check genome_pos is the expected value
    write_read_mapping(genome_pos, code_result_tab);
#endif

    code_result_idx = 0;
    while (code_result_tab[code_result_idx] != CODE_END) {
        int code_result = code_result_tab[code_result_idx];
        //LOG_DEBUG("code_result=%d\n", code_result);
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
                    //LOG_TRACE("set_variant early return\n");
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
    //LOG_TRACE("set_variant return\n");
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

#ifdef USE_INDEL
int get_read_update_positions(
        uint64_t * update_genome_position,
        dpu_result_out_t *result_tab, 
        int pos,
        genome_t *ref_genome,
        uint64_t genome_pos,
        int8_t *read,
        __attribute__((unused))int size_neighbour_in_symbols,
        bool * flag,
        bool debug,
        uint32_t * substCnt) {

    // run smith and waterman algorithm to find indels
    uint8_t code_result_tab[256];
    code_alignment(code_result_tab, result_tab[pos].score, &ref_genome->data[genome_pos], read, size_neighbour_in_symbols, flag);
    for(int read_pos = 0; read_pos < SIZE_READ; ++read_pos) {
        update_genome_position[read_pos] = 0;
    }
    if (code_result_tab[0] != CODE_ERR) {

        // array that contains for each read position, the genome position that it matches with
        // This is the genome position that will be updated in the frequency table
        // This genome position takes into account the shift due to possible indels found 
        // with smith-waterman algorithm
        int code_result_index = 0;
        int ref_pos = 0;
        int nbIndels = 0;
        bool ins = false;
        while (code_result_tab[code_result_index] != CODE_END) {
            int code_result = code_result_tab[code_result_index];
            int64_t pos_variant_read = code_result_tab[code_result_index + 1];
            /*printf("pos variant: %lu\n", pos_variant_read);*/
            int64_t pos_variant_genome = genome_pos + pos_variant_read;
            if (code_result == CODE_SUB) {
                // do nothing for substitution
                code_result_index += 3;
                (*substCnt)++;
                ref_pos++;
            }
            else if (code_result == CODE_INS) {
                ins = true;
                int64_t ps_var_genome = pos_variant_genome;
                int64_t ps_var_read = pos_variant_read;
                code_result_index += 2;

                while (code_result_tab[code_result_index] < 4) {
                    ps_var_read++;
                    code_result_index++;
                }

                while (ref_genome->data[ps_var_genome] == read[ps_var_read] && ps_var_genome 
                        && pos_variant_read) {
                    assert(ps_var_genome && ps_var_read && pos_variant_genome && pos_variant_read);
                    ps_var_genome--;
                    ps_var_read--;
                    pos_variant_genome--;
                    pos_variant_read--;
                }

                /*newvar->ref[ref_pos++] = nucleotide[ref_genome->data[pos_variant_genome] & 3];*/
                ref_pos++;

                // skip first value which should be the equivalent of first element in ref genome
                pos_variant_read++;
                while (pos_variant_read <= ps_var_read) {
                    // position should not be updated yet
                    if(update_genome_position[pos_variant_read] != 0) {
                        printf("Warning: duplicate update (Insertion) at position %lu. Current %lu\n", 
                                pos_variant_read, update_genome_position[pos_variant_read]);
                        fflush(stdout);
                        return -1;
                    }
                    update_genome_position[pos_variant_read++] = UINT64_MAX;
                }
                ++nbIndels;
            }
            else if (code_result == CODE_DEL) {

                int64_t ps_var_genome = pos_variant_genome;
                int64_t ps_var_read = pos_variant_read;
                code_result_index += 2;

                while (code_result_tab[code_result_index] < 4) {
                    ps_var_genome++;
                    code_result_index++;
                }

                while (ref_genome->data[ps_var_genome] == read[ps_var_read] && pos_variant_genome && ps_var_read) {
                    assert(ps_var_genome && ps_var_read && pos_variant_genome && pos_variant_read);
                    ps_var_read--;
                    ps_var_genome--;
                    pos_variant_genome--;
                    pos_variant_read--;
                }

                // on a deletion store the threshold to apply from the current read position
                assert(ps_var_genome > pos_variant_genome);
                if(pos_variant_read + 1 < SIZE_READ) {
                    // position should not be updated yet
                    if(update_genome_position[pos_variant_read+1] != 0) {
                        printf("Warning: duplicate update (Deletion) at position %lu. Current %lu\n", 
                                pos_variant_read+1, update_genome_position[pos_variant_read+1]);
                        fflush(stdout);
                        return -1;
                    }
                    /*assert(update_genome_position[pos_variant_read+1] == 0);*/
                    update_genome_position[pos_variant_read + 1] = ps_var_genome - pos_variant_genome;
                }
                while (pos_variant_genome <= ps_var_genome) {
                    pos_variant_genome++;
                    ref_pos++;
                }
                pos_variant_genome -= ref_pos;
                ++nbIndels;
            }
            else
                assert(0);
        }

        // debug prints
        if(nbIndels && debug)
            printf("SW algorithm (nbIndels %d) ins %d:\n", nbIndels, ins);
        int64_t curr_pos = genome_pos;
        for(int read_pos = 0; read_pos < SIZE_READ; ++read_pos) {
            switch(update_genome_position[read_pos]) {
                case 0:
                    update_genome_position[read_pos] = curr_pos++;
                    if(nbIndels && debug)
                        printf(" ");
                    break;
                case UINT64_MAX:
                    if(nbIndels && debug)
                        printf("I");
                    break;
                default:
                    if(nbIndels && debug) {
                        for(uint64_t print_index = 0; print_index < update_genome_position[read_pos]; ++print_index)
                            printf("D");
                    }
                    curr_pos += update_genome_position[read_pos];
                    update_genome_position[read_pos] = curr_pos++;
            }
        }

        if(nbIndels && debug) {
            printf("\n");
            fflush(stdout);
        }
        return nbIndels;
    }
    else
        assert(0);

    return false;
}
#endif


static pthread_mutex_t freq_table_mutex;

/**
 * function to update frequency table used for variant calling
 **/
bool update_frequency_table(
        genome_t *ref_genome, 
        dpu_result_out_t *result_tab, 
        int8_t *reads_buffer, 
        float *reads_quality_buffer, 
        int pos,
        float mapq,
        __attribute__((unused))int size_neighbour_in_symbols) {

    struct frequency_info **frequency_table = get_frequency_table();
    uint64_t genome_pos = ref_genome->pt_seq[result_tab[pos].coord.seq_nr] + result_tab[pos].coord.seed_nr; 
    int num = result_tab[pos].num;
    int8_t *read = reads_buffer + (num * SIZE_READ);
    float *read_quality = reads_quality_buffer + (num/2 * SIZE_READ);
    //TODO: read the quality in the correct order (inverted or not)
    bool inv = num & 1;
    //TODO: assume no offset here

#ifdef USE_INDEL

    static bool debug = false; 
    static char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    uint64_t update_genome_position[SIZE_READ];
    uint32_t substCnt = 0;
    flag_dbg = false;

    // for simplicity put all this in a critical section protected by a mutex
    // since the frequency table is shared (but inefficient)
    
    pthread_mutex_lock(&freq_table_mutex);
    int nbIndels = get_read_update_positions(update_genome_position, result_tab, pos,
            ref_genome, genome_pos, read, size_neighbour_in_symbols, &flag_dbg, debug, &substCnt);
    bool hasIndel = nbIndels > 0;

    // debug prints
    if(hasIndel && debug) {

        printf("Read:\n");
        for(int k = 0; k < SIZE_READ; ++k) {
            printf("%c", nucleotide[read[k]]);
        }
        printf("\ngenome (pos %u:%u):\n", result_tab[pos].coord.seed_nr, result_tab[pos].coord.seed_nr+SIZE_READ);
        for(uint64_t k = genome_pos; k < genome_pos + SIZE_READ; ++k) {
            printf("%c", nucleotide[ref_genome->data[k]]);
        }
        printf("\nupdate pos:\n");
        uint64_t lastpos = 0;
        for(uint64_t k = 0; k < SIZE_READ; ++k) {
            if(k && update_genome_position[k] != lastpos+1) {
                if(update_genome_position[k] == UINT64_MAX)
                    printf("X");
                /*printf("No update at position %lu\n", k);*/
                else if(lastpos == UINT64_MAX)
                    /*printf("New start at pos %lu = %lu / %c\n", k, update_genome_position[k], nucleotide[ref_genome->data[update_genome_position[k]]]);*/
                    printf("%c", nucleotide[ref_genome->data[update_genome_position[k]]]);
                else
                    /*printf("Change at pos %lu, diff %ld, %c\n", k, update_genome_position[k] - lastpos, nucleotide[ref_genome->data[update_genome_position[k]]]);*/
                    printf("%c", nucleotide[ref_genome->data[update_genome_position[k]]]);
            }
            /*else if(nucleotide[ref_genome->data[update_genome_position[k]]] != nucleotide[read[k]]) {*/
            /*printf(RED "%c" RESET, nucleotide[ref_genome->data[update_genome_position[k]]]);*/
            /*}*/
            else
                printf("%c", nucleotide[ref_genome->data[update_genome_position[k]]]);
            lastpos = update_genome_position[k];
        }

        printf("\nsubst:\n");
        for(uint64_t k = 0; k < SIZE_READ; ++k) {
            if(update_genome_position[k] == UINT64_MAX) {
                printf(" ");
                continue;
            }
            else if(nucleotide[ref_genome->data[update_genome_position[k]]] != nucleotide[read[k]]) {
                printf("U");
                substCnt++;
            }
            else
                printf(" ");
        }
        printf("\n\n");
        fflush(stdout);
        assert(!result_tab[pos].coord.nodp);
    }
    else if(debug) {
        if(!result_tab[pos].coord.nodp) {
            printf("\nWarning: odpd result with no indels detected (flag = %d, subst cnt %u)):\n", flag_dbg, substCnt);
            printf("Read:\n");
            for(int k = 0; k < SIZE_READ; ++k) {
                printf("%c", nucleotide[read[k]]);
            }
            printf("\ngenome (pos %u:%u):\n", result_tab[pos].coord.seed_nr, result_tab[pos].coord.seed_nr+SIZE_READ);
            for(uint64_t k = genome_pos; k < genome_pos + SIZE_READ; ++k) {
                printf("%c", nucleotide[ref_genome->data[k]]);
            }
            printf("\n\n");
            fflush(stdout);
        }
    }

    /*pthread_mutex_lock(&freq_table_mutex);*/
    // for the moment support only one indel otherwise we have some issue FIXME
    if(substCnt <= MAX_SUBSTITUTION && (hasIndel || result_tab[pos].coord.nodp) && nbIndels >= 0) {
        for(uint64_t k = 0; k < SIZE_READ; ++k) {
            uint64_t update_genome_pos = update_genome_position[k];
            if(update_genome_pos < genome_get()->fasta_file_size) {
                frequency_table[read[k]][update_genome_pos].freq += mapq * read_quality[inv ? SIZE_READ - k - 1 : k];
                /*frequency_table[read[j]][genome_pos+j].score += result_tab[pos].score;*/
                frequency_table[read[k]][update_genome_pos].score++;
            }
            else if (update_genome_pos != UINT64_MAX)
                printf("WARNING: genome update position computed is wrong %lu\n", update_genome_pos);
        }
    }
    /*fflush(stdout);*/
    pthread_mutex_unlock(&freq_table_mutex);
    return hasIndel;

#else
    pthread_mutex_lock(&freq_table_mutex);
    for(int j = 0; j < SIZE_READ; ++j) {
        if(genome_pos + j < genome_get()->fasta_file_size) {
            frequency_table[read[j]][genome_pos+j].freq += mapq * read_quality[inv ? SIZE_READ - j - 1 : j];
            frequency_table[read[j]][genome_pos+j].score++;
        }
        else
            printf("WARNING: reads matched at position that exceeds genome size\n");
    }
    pthread_mutex_unlock(&freq_table_mutex);
    return false;
#endif
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
    float *reads_quality_buffer;
    genome_t *ref_genome;
    FILE *fpe1;
    FILE *fpe2;
} process_read_arg_t;

static uint64_t nr_reads_total = 0ULL;
static uint64_t nr_reads_total_from_dpus = 0ULL;
static uint64_t nr_reads_non_mapped = 0ULL;
static uint64_t nr_reads_with_indels = 0ULL;
static pthread_mutex_t nr_reads_mutex = PTHREAD_MUTEX_INITIALIZER;

#define MISMATCH_COUNT(X) (X.score / 10) 
#define INVALID_SCORE 1000

static void keep_best_2_scores(unsigned score, unsigned* P1, unsigned *P2, unsigned x1, unsigned x2, unsigned* best_score) {

  if(score < best_score[0]) {

    // move current to next position
    best_score[1] = best_score[0];
    P1[1] = P1[0];
    P2[1] = P2[0];
    // update first position
    best_score[0] = score;
    P1[0] = x1;
    P2[0] = x2;
  }
  else if (score < best_score[1]) {

    // update second position
    best_score[1] = score;
    P1[1] = x1;
    P2[1] = x2;
  }
}

static unsigned get_nb_scores(unsigned int * best_score) {

    unsigned np = 0;
    if(best_score[0] < INVALID_SCORE) {
        np++;
        if(best_score[1] < INVALID_SCORE) np++;
    }
    return np;
}

/*#define USE_MAPQ_SCORE*/
static void do_process_read(process_read_arg_t *arg)
{
    const unsigned int nb_match = arg->nb_match;
    dpu_result_out_t *result_tab = arg->result_tab;
    int round = arg->round;
    int8_t *reads_buffer = arg->reads_buffer;
    float *reads_quality_buffer = arg->reads_quality_buffer;
    genome_t *ref_genome = arg->ref_genome;
    FILE *fpe1 = arg->fpe1;
    FILE *fpe2 = arg->fpe2;
    unsigned int size_neighbour_in_symbols = (SIZE_NEIGHBOUR_IN_BYTES - DELTA_NEIGHBOUR(round)) * 4;
    /*printf("size_neighbour_in_symbols : %u", size_neighbour_in_symbols);*/

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
      unsigned int P1[2];
      unsigned int P2[2];
      unsigned int pos1, pos2, t1, t2;
      unsigned int best_score[2] = { 1000, 1000 };
      /*unsigned int best_score_all = 1000;*/
      // test all significant pairs of reads (0,3) & (1,2)
      for (unsigned int x1 = i; x1 < j; x1++) {
        t1 = result_tab[x1].num % 4;
        pos1 = result_tab[x1].coord.seed_nr;
        for (unsigned int x2 = i + 1; x2 < j; x2++) {
          pos2 = result_tab[x2].coord.seed_nr;
          t2 = result_tab[x2].num % 4;
          if (t1 + t2 == 3) // select significant pair
          {
            if ((abs((int)pos2 - (int)pos1) > READ_DIST_LOWER_BOUND && (abs((int)pos2 - (int)pos1) < READ_DIST_UPPER_BOUND))) {
              // update if this is one of the two best scores
              keep_best_2_scores(result_tab[x1].score + result_tab[x2].score, P1, P2, x1, x2, best_score);
            }
          }
        }
      }

      bool update = false;
      bool hasIndel = false;

      unsigned np = get_nb_scores(best_score);
      if (np > 0) {

        if(np == 2) {

          // found at least 2 matching pairs of positions. Check the delta between the two pairs to 
          // decide whether we should keep the best pair
          int delta = abs((int)(MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]])) 
              - (int)(MISMATCH_COUNT(result_tab[P1[1]]) + MISMATCH_COUNT(result_tab[P2[1]])));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) 
            + MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((2 * (MAX_SUBSTITUTION + 1)) - delta);
          if(delta_corrected < 0) {
            printf("WARNING: negative delta for square root %d\n", delta_corrected);
          }
          else if(delta > DIST_PAIR_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta > DIST_PAIR_THRESHOLD) {
#endif
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, size_neighbour_in_symbols);
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }
        else if(np) { // only one result, take it
          int delta = abs((int)(MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]])) - (2 * (MAX_SUBSTITUTION + 1)));
          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((2 * (MAX_SUBSTITUTION + 1)) - delta);
          if(delta_corrected < 0) {
            printf("WARNING: negative delta (np == 1) for square root %d\n", delta_corrected);
          }
          else if(delta > DIST_PAIR_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta > DIST_PAIR_THRESHOLD) {
#endif
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, size_neighbour_in_symbols);
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }
      }
      if(true) {

        // check mapping of R1 and R2 independently
        unsigned int best_score_R1[2] = { 1000, 1000 };
        unsigned int best_score_R2[2] = { 1000, 1000 };
        P1[0] = 0;
        P2[0] = 0;
        P1[1] = 0;
        P2[1] = 0;
        for (unsigned int read = i; read < j; read++) {
          unsigned t1 = result_tab[read].num % 4;
          if(t1 < 2) { // PE1 or RPE1
            keep_best_2_scores(result_tab[read].score, P1, P2, read, 0, best_score_R1);
          }
          else { // PE2 or RPE2
            keep_best_2_scores(result_tab[read].score, P1, P2, 0, read, best_score_R2);
          }
        }

        unsigned np1 = get_nb_scores(best_score_R1), np2 = get_nb_scores(best_score_R2);
        if(np1 == 2) {

          int delta = abs((int)MISMATCH_COUNT(result_tab[P1[0]]) - (int)MISMATCH_COUNT(result_tab[P1[1]]));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            printf("WARNING: negative delta (np1 == 2) for square root %d\n", delta_corrected);
          }
          else if(delta > DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta > DIST_SINGLE_THRESHOLD) {
#endif
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }
        else if(np1) {
          int delta = abs((int)MISMATCH_COUNT(result_tab[P1[0]]) - (MAX_SUBSTITUTION + 1));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            printf("WARNING: negative delta (np1 == 1) for square root %d\n", delta_corrected);
          }
          else if(delta > DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta > DIST_SINGLE_THRESHOLD) {
#endif
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }

        if(np2 == 2) {

          int delta = abs((int)MISMATCH_COUNT(result_tab[P2[0]]) - (int)MISMATCH_COUNT(result_tab[P2[1]]));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            printf("WARNING: negative delta (np2 == 2) for square root %d, %d %d %d %d\n", 
                delta_corrected, MISMATCH_COUNT(result_tab[P2[0]]), MISMATCH_COUNT(result_tab[P2[1]]), MAX_SUBSTITUTION + 1, delta);
          }
          else if(delta > DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta > DIST_SINGLE_THRESHOLD) {
#endif

            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }
        else if(np2) {
          int delta = abs((int)MISMATCH_COUNT(result_tab[P2[0]]) - (MAX_SUBSTITUTION + 1));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            printf("WARNING: negative delta (np2 == 1) for square root %d\n", delta_corrected);
          }
          else if (delta > DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if (delta > DIST_SINGLE_THRESHOLD) {
#endif
            hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, size_neighbour_in_symbols);
            update = true;
          }
        }
        if(!update) {
            pthread_mutex_lock(&nr_reads_mutex);
            nr_reads_non_mapped++;
            pthread_mutex_unlock(&nr_reads_mutex);
            add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
        }
        if(hasIndel)
            nr_reads_with_indels++;
        pthread_mutex_lock(&nr_reads_mutex);
        nr_reads_total_from_dpus++;
        pthread_mutex_unlock(&nr_reads_mutex);
      }
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
    float *reads_quality_buffer = get_reads_quality_buffer(pass_id);
    acc_results_t acc_res = accumulate_get_result(pass_id);
    nr_reads_total += get_reads_in_buffer(pass_id) / 4;

    curr_match = 0;

    args.nb_match = acc_res.nb_res;
    args.result_tab = acc_res.results;
    args.round = round;
    args.reads_buffer = reads_buffer;
    args.reads_quality_buffer = reads_quality_buffer;
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
#if DEBUG_READ_MAPPING
    open_sam_file();
#endif
    genome_t *ref_genome = genome_get();
    args.ref_genome = ref_genome;

    assert(pthread_mutex_init(&curr_match_mutex, NULL) == 0);
    assert(pthread_mutex_init(&non_mapped_mutex, NULL) == 0);
    assert(pthread_mutex_init(&freq_table_mutex, NULL) == 0);
    assert(pthread_barrier_init(&barrier, NULL, PROCESS_READ_THREAD) == 0);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_create(&thread_id[each_thread], NULL, process_read_thread_fct, &args) == 0);
    }
}

void process_read_free()
{
#if DEBUG_READ_MAPPING
    close_sam_file();
#endif
    stop_threads = true;
    pthread_barrier_wait(&barrier);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_join(thread_id[each_thread], NULL) == 0);
    }

    assert(pthread_barrier_destroy(&barrier) == 0);
    assert(pthread_mutex_destroy(&curr_match_mutex) == 0);
    assert(pthread_mutex_destroy(&non_mapped_mutex) == 0);
    assert(pthread_mutex_destroy(&freq_table_mutex) == 0);
    fflush(stdout);
    fprintf(stderr, "%% reads non mapped: %f%%\n", (float)nr_reads_non_mapped * 100.0 / (float)nr_reads_total_from_dpus);
    fprintf(stderr, "%% reads with indels: %f%%\n", (float)nr_reads_with_indels * 100.0 / (float)(nr_reads_total_from_dpus - nr_reads_non_mapped));
    fprintf(stderr, "%% Total reads from dpus: %ld%%\n", nr_reads_total_from_dpus);
    fprintf(stderr, "%% Total reads: %ld%%\n", nr_reads_total);
}
