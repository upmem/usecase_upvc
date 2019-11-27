/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "compare.h"
#include "genome.h"
#include "parse_args.h"
#include "upvc.h"
#include "upvc_dpu.h"
#include "vartree.h"

#include "common.h"

#define SIZE_INSERT_MEAN (400)
#define SIZE_INSERT_STD (3 * 50)

static int cmpresult(const void *a, const void *b)
{
    dpu_result_out_t *A = (dpu_result_out_t *)a;
    dpu_result_out_t *B = (dpu_result_out_t *)b;
    if (A->num == B->num) {
        if (A->score > B->score) {
            return 1;
        } else {
            return -1;
        }
    }
    if (A->num > B->num) {
        return 1;
    } else {
        return -1;
    }
}

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

static void set_variant(dpu_result_out_t result_match, genome_t *ref_genome, int8_t *reads_buffer, variant_tree_t **variant_list,
                        int *substitution_list, int8_t *mapping_coverage, unsigned int size_neighbour_in_symbols)
{
    uint32_t code_result_idx;
    uint8_t code_result_tab[256];
    int8_t *read;
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    uint64_t genome_pos = ref_genome->pt_seq[result_match.coord.seq_nr] + result_match.coord.seed_nr;
    int size_read = SIZE_READ;

    /* Update "mapping_coverage" with the number of reads that match at this position of the genome */
    for (int i = 0; i < size_read; i++) {
        mapping_coverage[genome_pos + i] += 1;
    }

    /* Get the differences betweend the read and the sequence of the reference genome that match */
    read = &reads_buffer[result_match.num * size_read];
    code_alignment(code_result_tab, result_match.score, &ref_genome->data[genome_pos], read, size_neighbour_in_symbols);

    code_result_idx = 0;
    while (code_result_tab[code_result_idx] != CODE_END) {
        int code_result = code_result_tab[code_result_idx];
        int64_t pos_variant_read = code_result_tab[code_result_idx + 1];
        int64_t pos_variant_genome = genome_pos + pos_variant_read;
        if (code_result == CODE_SUB) {
            /* SNP = 0,1,2,3  (code A,C,T,G) */
            int snp = code_result_tab[code_result_idx + 2];

            /* substitution_list[pos_variant_genome] =  32 bits integer containning 4 counter of 8 bits */
            substitution_list[pos_variant_genome] += 1 << (snp * 8);

            code_result_idx += 3;
        } else if (code_result == CODE_INS) {
            int64_t ps_var_genome = pos_variant_genome;
            int64_t ps_var_read = pos_variant_read;
            variant_t *newvar = (variant_t *)malloc(sizeof(variant_t));

            newvar->offset = ref_genome->pt_seq[result_match.coord.seq_nr];
            newvar->chr = ref_genome->seq_name[result_match.coord.seq_nr];
            newvar->depth = 1;
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

            newvar->pos = pos_variant_genome;
            newvar->ref[0] = nucleotide[ref_genome->data[pos_variant_genome] & 3];
            newvar->ref[1] = '\0';

            {
                int k = 0;
                while (pos_variant_read <= ps_var_read) {
                    newvar->alt[k] = nucleotide[read[pos_variant_read] & 3];
                    k++;
                    pos_variant_read++;
                }

                newvar->alt[k] = '\0';
                insert_variants(variant_list, newvar);
            }
        } else if (code_result == CODE_DEL) {
            int64_t ps_var_genome = pos_variant_genome;
            int64_t ps_var_read = pos_variant_read;
            variant_t *newvar = (variant_t *)malloc(sizeof(variant_t));

            newvar->offset = ref_genome->pt_seq[result_match.coord.seq_nr];
            newvar->chr = ref_genome->seq_name[result_match.coord.seq_nr];
            newvar->depth = 1;
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

            newvar->pos = pos_variant_genome;
            newvar->alt[0] = nucleotide[ref_genome->data[pos_variant_genome] & 3];
            newvar->alt[1] = '\0';

            {
                int k = 0;
                while (pos_variant_genome <= ps_var_genome) {
                    newvar->ref[k] = nucleotide[ref_genome->data[pos_variant_genome] & 3];
                    k++;
                    pos_variant_genome++;
                }

                newvar->ref[k] = '\0';
                insert_variants(variant_list, newvar);
            }
        }
    }
}

static void add_to_non_mapped_read(int numread, int round, FILE *fpe1, FILE *fpe2, int8_t *reads_buffer)
{
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
}

static bool compute_read_pair(int type1, int type2, int *offset, int *nbread, dpu_result_out_t *result_tab, int *result_tab_tmp,
    int numpair, int round, FILE *fpe1, FILE *fpe2, int8_t *reads_buffer, bool *result_found)
{
    for (int pa = offset[type1]; pa < offset[type1] + nbread[type1]; pa++) {
        for (int pb = offset[type2]; pb < offset[type2] + nbread[type2]; pb++) {
            if (check_pair(result_tab[pa], result_tab[pb])) {
                if (*result_found) {
                    add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
                    return false;
                }
                result_tab_tmp[0] = pa;
                result_tab_tmp[1] = pb;
                *result_found = true;
            }
        }
    }
    return true;
}

int accumulate_read(
    dpu_result_out_t **result_tab, unsigned int *result_tab_nb_read, unsigned int dpu_offset, times_ctx_t *times_ctx)
{
    double t1, t2;
    unsigned int nb_match = *result_tab_nb_read;
    unsigned int nb_total_read = 0;

    t1 = my_clock();

    for (unsigned int numdpu = dpu_offset; numdpu < dpu_offset + get_nb_dpus_per_run(); numdpu++) {
        unsigned int k = 0;
        int num_read;
        do {
            num_read = read_out_num(numdpu, k);
            k++;
        } while (num_read != -1);
        nb_total_read += (k - 1);
    }
    *result_tab_nb_read = *result_tab_nb_read + nb_total_read;
    *result_tab = realloc(*result_tab, *result_tab_nb_read * sizeof(dpu_result_out_t));

    for (unsigned int numdpu = dpu_offset; numdpu < dpu_offset + get_nb_dpus_per_run(); numdpu++) {
        int k = 0;
        int num_read = read_out_num(numdpu, k);

        while (num_read != -1) {
            dpu_result_coord_t coord = read_out_coord(numdpu, k);

            (*result_tab)[nb_match].num = num_read;
            (*result_tab)[nb_match].coord.seq_nr = coord.seq_nr;
            (*result_tab)[nb_match].coord.seed_nr = coord.seed_nr;
            (*result_tab)[nb_match].score = read_out_score(numdpu, k);
            nb_match++;

            num_read = read_out_num(numdpu, ++k);
        }
    }
    assert(nb_match == *result_tab_nb_read);
    qsort(*result_tab, nb_match, sizeof(dpu_result_out_t), cmpresult);

    t2 = my_clock();
    times_ctx->acc_read = t2 - t1;
    times_ctx->tot_acc_read += t2 - t1;

    return (int)nb_total_read;
}

int process_read(genome_t *ref_genome, int8_t *reads_buffer, variant_tree_t **variant_list, int *substitution_list,
    int8_t *mapping_coverage, dpu_result_out_t *result_tab, unsigned int result_tab_nb_read, FILE *fpe1, FILE *fpe2, int round,
    times_ctx_t *times_ctx)
{
    double t1, t2;
    int nb_match = result_tab_nb_read;
    int nb_read_map = 0;
    int offset[4]; /* offset in result_tab of the first read  */
    int nbread[4]; /* number of reads                        */
    unsigned int score[4];
    int result_tab_tmp[2];
    unsigned int size_neighbour_in_symbols = (SIZE_NEIGHBOUR_IN_BYTES - DELTA_NEIGHBOUR(round)) * 4;

    t1 = my_clock();

    /*
     * The number of a pair is given by "num_read / 4 " (see dispatch_read function)
     * Their type is given by their offset (see dispatch_read function)
     * type = num_read%4 == 0 ==> read 1
     *                      1 ==> read 1 complement
     *                      2 ==> read 2
     *                      3 ==> read 2 complement
     * The read pair to consider are [0, 3] and [1, 2].
     */

    int i = 0;
    while (i < nb_match) {
        int numpair = result_tab[i].num / 4;
        bool result_found = false;
        for (int k = 0; k < 4; k++) {
            score[k] = 1000;
            offset[k] = -1;
            nbread[k] = 0;
        }
        while ((i < nb_match) && (numpair == result_tab[i].num / 4)) {
            int type = result_tab[i].num % 4;
            if (result_tab[i].score < score[type]) {
                score[type] = result_tab[i].score;
                offset[type] = i;
                nbread[type] = 1;
            } else if (result_tab[i].score == score[type]) {
                nbread[type]++;
            }
            i++;
        }

        /* Compute a [0, 3] read pair */
        if (!compute_read_pair(0, 3, offset, nbread, result_tab, result_tab_tmp, numpair, round, fpe1, fpe2, reads_buffer,
                &result_found)) {
            continue;
        }

        /* Compute a [1, 2] read pair */
        if (!compute_read_pair(2, 1, offset, nbread, result_tab, result_tab_tmp, numpair, round, fpe1, fpe2, reads_buffer,
                &result_found)) {
            continue;
        }

        if (result_found) { /* Mapped reads*/
            set_variant(result_tab[result_tab_tmp[0]], ref_genome, reads_buffer, variant_list, substitution_list,
                mapping_coverage, size_neighbour_in_symbols);
            set_variant(result_tab[result_tab_tmp[1]], ref_genome, reads_buffer, variant_list, substitution_list,
                mapping_coverage, size_neighbour_in_symbols);
            nb_read_map += 2;
        } else {
            add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
        }
    }

    free(result_tab);

    t2 = my_clock();
    times_ctx->process_read = t2 - t1;
    times_ctx->tot_process_read += t2 - t1;

    return nb_read_map;
}
