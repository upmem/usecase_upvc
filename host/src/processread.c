/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "accumulateread.h"
#include "compare.h"
#include "genome.h"
#include "getread.h"
#include "upvc.h"
#include "vartree.h"

#define SIZE_INSERT_MEAN (400)
#define SIZE_INSERT_STD (3 * 50)

static int check_pair(dpu_result_out_t A1, dpu_result_out_t A2)
{
    int pos1, pos2;
    int size_insert_min = SIZE_INSERT_MEAN - SIZE_INSERT_STD;
    int size_insert_max = SIZE_INSERT_MEAN + SIZE_INSERT_STD;
    int size_read = SIZE_READ;

    if (A1.coord.seq_nr == A2.coord.seq_nr) {
        pos1 = A1.coord.seed_nr;
        pos2 = A2.coord.seed_nr;

        /* Check the distance between the 2 reads */
        if ((abs(pos2 - pos1 + size_read) > size_insert_min) && (abs(pos2 - pos1 + size_read) < size_insert_max)) {
            return 1;
        }
    }
    return 0;
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

static void set_variant(dpu_result_out_t result_match, genome_t *ref_genome, int8_t *reads_buffer, unsigned int size_neighbour_in_symbols)
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
        newvar->depth=1;
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

static pthread_mutex_t non_mapped_mutex;
#if NB_ROUND > 1
static void add_to_non_mapped_read(int numread, int round, FILE *fpe1, FILE *fpe2, int8_t *reads_buffer)
{
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
#else
static void add_to_non_mapped_read(__attribute__((unused)) int numread, __attribute__((unused)) int round,
    __attribute__((unused)) FILE *fpe1, __attribute__((unused)) FILE *fpe2, __attribute__((unused)) int8_t *reads_buffer)
{
}
#endif

static bool compute_read_pair(
    int type1, int type2, int *offset, int *nbread, dpu_result_out_t *result_tab, int *pa_found, int *pb_found, bool *found)
{
    for (int pa = offset[type1]; pa < offset[type1] + nbread[type1]; pa++) {
        for (int pb = offset[type2]; pb < offset[type2] + nbread[type2]; pb++) {
            if (check_pair(result_tab[pa], result_tab[pb])) {
                if (*found)
                    return false;
                *found = true;
                *pa_found = pa;
                *pb_found = pb;
            }
        }
    }
    return true;
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

#define PROCESS_READ_THREAD (8)
typedef struct {
    const unsigned int nb_match;
    dpu_result_out_t *result_tab;
    int round;
    int8_t *reads_buffer;
    genome_t *ref_genome;
    FILE *fpe1;
    FILE *fpe2;
} process_read_arg_t;

void * process_read_thread_fct(void *arg) {
    const unsigned int nb_match = ((process_read_arg_t *)arg)->nb_match;
    dpu_result_out_t *result_tab = ((process_read_arg_t *)arg)->result_tab;
    int round = ((process_read_arg_t *)arg)->round;
    int8_t *reads_buffer = ((process_read_arg_t *)arg)->reads_buffer;
    genome_t *ref_genome = ((process_read_arg_t *)arg)->ref_genome;
    FILE * fpe1 = ((process_read_arg_t *)arg)->fpe1;
    FILE * fpe2 = ((process_read_arg_t *)arg)->fpe2;
    unsigned int size_neighbour_in_symbols = (SIZE_NEIGHBOUR_IN_BYTES - DELTA_NEIGHBOUR(round)) * 4;

    /*
     * The number of a pair is given by "num_read / 4 " (see dispatch_read function)
     * Their type is given by their offset (see dispatch_read function)
     * type = num_read%4 == 0 ==> read 1
     *                      1 ==> read 1 complement
     *                      2 ==> read 2
     *                      3 ==> read 2 complement
     * The read pair to consider are [0, 3] and [1, 2].
     */

    while (true) {
        int offset[4] = { -1, -1, -1, -1 }; /* offset in result_tab of the first read  */
        int nbread[4] = { 0, 0, 0, 0 }; /* number of reads                        */
        unsigned int score[4] = { 1000, 1000, 1000, 1000 };

        unsigned int i;
        if ((i = acquire_curr_match()) >= nb_match) {
            release_curr_match(i);
            return NULL;
        }
        int numpair = result_tab[i].num / 4;
        unsigned int j = i;
        while ((j < nb_match) && (numpair == result_tab[j].num / 4)) {
            j++;
        }
        release_curr_match(j);

        for (; i < j; i++) {
            int type = result_tab[i].num % 4;
            if (result_tab[i].score < score[type]) {
                score[type] = result_tab[i].score;
                offset[type] = i;
                nbread[type] = 1;
            } else if (result_tab[i].score == score[type]) {
                nbread[type]++;
            }
        }

        bool found = false;
        int pa;
        int pb;
        /* Compute a [0, 3] read pair */
        if (!compute_read_pair(0, 3, offset, nbread, result_tab, &pa, &pb, &found)) {
            goto add_to_non_mapped_read_label;
        }
        /* Compute a [1, 2] read pair */
        if (!compute_read_pair(2, 1, offset, nbread, result_tab, &pa, &pb, &found)) {
            goto add_to_non_mapped_read_label;
        }

        if (found) {
            set_variant(result_tab[pa], ref_genome, reads_buffer, size_neighbour_in_symbols);
            set_variant(result_tab[pb], ref_genome, reads_buffer, size_neighbour_in_symbols);
            continue;
        }

    add_to_non_mapped_read_label:
        add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
    }
}

void process_read(FILE *fpe1, FILE *fpe2, int round, unsigned int pass_id)
{
    unsigned int nb_match;
    dpu_result_out_t *result_tab;
    int8_t *reads_buffer = get_reads_buffer(pass_id);
    genome_t *ref_genome = genome_get();

    accumulate_get_result(pass_id, &nb_match, &result_tab);

    pthread_mutex_init(&curr_match_mutex, NULL);
    pthread_mutex_init(&non_mapped_mutex, NULL);
    curr_match = 0;

    process_read_arg_t args = {
        .nb_match = nb_match,
        .result_tab = result_tab,
        .round = round,
        .reads_buffer = reads_buffer,
        .ref_genome = ref_genome,
        .fpe1 = fpe1,
        .fpe2 = fpe2,
    };
    int ret;
    pthread_t thread_id[PROCESS_READ_THREAD];
    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD; each_thread++) {
        ret = pthread_create(&thread_id[each_thread], NULL, process_read_thread_fct, &args);
        assert(ret == 0);
    }
    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD; each_thread++) {
        ret = pthread_join(thread_id[each_thread], NULL);
        assert(ret == 0);
    }

    pthread_mutex_destroy(&curr_match_mutex);
    pthread_mutex_destroy(&non_mapped_mutex);
    free(result_tab);
}
