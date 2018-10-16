#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "upvc_dpu.h"
#include "compare.h"
#include "vartree.h"
#include "genome.h"
#include "upvc.h"

// alignement
typedef struct {
        int num_read;
        int coord_seq;
        int coord_pos;
        int score;
} align_t;

#define  SIZE_INSERT_MEAN (400)
#define  SIZE_INSERT_STD (3 * 50)

static int cmpalign(const void * a, const void * b)
{
        align_t *A = (align_t *) a;
        align_t *B = (align_t *) b;
        if (A->num_read == B->num_read) {
                if (A->score > B->score) {
                        return 1;
                } else {
                        return -1;
                }
        }
        if (A->num_read > B->num_read) {
                return 1;
        } else {
                return -1;
        }
}

static int check_pair(align_t A1, align_t A2, reads_info_t *reads_info)
{
        int pos1, pos2;
        int size_insert_min = SIZE_INSERT_MEAN - SIZE_INSERT_STD;
        int size_insert_max = SIZE_INSERT_MEAN + SIZE_INSERT_STD;
        int size_read = reads_info->size_read;

        if (A1.coord_seq == A2.coord_seq) {
                pos1 = A1.coord_pos;
                pos2 = A2.coord_pos;

                /* Check the distance between the 2 reads */
                if ( (abs(pos2 - pos1 + size_read) > size_insert_min)
                     &&
                     (abs(pos2 - pos1 + size_read) < size_insert_max) ) {
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

static int code_alignment(int8_t *code, int score, int8_t *gen, int8_t *read, reads_info_t *reads_info)
{
        int code_idx, computed_score, backtrack_idx;
        int size_read = reads_info->size_read;
        int size_neighbour = reads_info->size_neighbour_in_32bits_words;
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
        backtrack_idx = DPD(&gen[SIZE_SEED], &read[SIZE_SEED], backtrak, reads_info);

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
                                while ( (backtrak[backtrack_idx].type == CODE_DEL)
                                        &&
                                        (backtrack_jx == backtrak[backtrack_idx].jx) ) {
                                        code[code_idx++] = gen[backtrak[backtrack_idx].ix + SIZE_SEED] & 3;
                                        backtrack_idx--;
                                }
                        } else {
                                int backtrack_ix = backtrak[backtrack_idx].ix;
                                code[code_idx++] = CODE_INS;
                                code[code_idx++] = backtrak[backtrack_idx].jx + SIZE_SEED - 1;
                                code[code_idx++] = read[backtrak[backtrack_idx].jx + SIZE_SEED];
                                backtrack_idx--;
                                while ( (backtrak[backtrack_idx].type == CODE_INS)
                                        &&
                                        (backtrack_ix == backtrak[backtrack_idx].ix) ) {
                                        code[code_idx++] = read[backtrak[backtrack_idx].jx + SIZE_SEED];
                                        backtrack_idx--;
                                }
                        }
                }
        }
        code[code_idx++] = CODE_END;
        return code_idx;
}

static void set_variant(align_t A,
                        genome_t *ref_genome,
                        int8_t *reads_buffer,
                        variant_tree_t **variant_list,
                        int *substitution_list,
                        int8_t *mapping_coverage,
                        reads_info_t *reads_info)
{
        int code_align_idx;
        int8_t code_align_tab[256];
        int8_t *read;
        char nucleotide[4] = {'A','C','T','G'};
        int genome_pos = ref_genome->pt_seq[A.coord_seq] + A.coord_pos;
        int size_read = reads_info->size_read;

        /* Update "mapping_coverage" withe the number of reads that match at this position of the genome */
        for (int i = 0; i < size_read; i++) {
                mapping_coverage[genome_pos + i] += 1;
        }

        /* Get the differences betweend the read and the sequence of the reference genome that match */
        read = &reads_buffer[A.num_read * size_read];
        code_alignment(code_align_tab, A.score, &ref_genome->data[genome_pos], read, reads_info);

        code_align_idx = 0;
        while (code_align_tab[code_align_idx] != CODE_END) {
                int code_align = code_align_tab[code_align_idx];
                int pos_variant_read = code_align_tab[code_align_idx + 1];
                int pos_variant_genome = genome_pos + pos_variant_read;
                if (code_align == CODE_SUB) {
                        /* SNP = 0,1,2,3  (code A,C,T,G) */
                        int snp = code_align_tab[code_align_idx + 2];

                        /* substitution_list[pos_variant_genome] =  32 bits integer containning 4 counter of 8 bits */
                        substitution_list[pos_variant_genome] += 1 << (snp * 8);

                        code_align_idx += 3;
                } else if (code_align == CODE_INS) {
                        int ps_var_genome = pos_variant_genome;
                        int ps_var_read = pos_variant_read;
                        variant_t *newvar = (variant_t*) malloc(sizeof(variant_t));

                        newvar->offset = ref_genome->pt_seq[A.coord_seq];
                        newvar->chr = ref_genome->seq_name[A.coord_seq];
                        newvar->depth = 1;
                        code_align_idx += 2;

                        while (code_align_tab[code_align_idx] < 4) {
                                ps_var_read++;
                                code_align_idx++;
                        }

                        while (ref_genome->data[ps_var_genome] == read[ps_var_read]) {
                                ps_var_genome--;
                                ps_var_read--;
                                pos_variant_genome--;
                                pos_variant_read--;
                        }

                        newvar->pos = pos_variant_genome;
                        newvar->ref[0] = nucleotide[ref_genome->data[pos_variant_genome]&3];
                        newvar->ref[1] = '\0';

                        {
                                int k = 0;
                                while (pos_variant_read <= ps_var_read) {
                                        newvar->alt[k] = nucleotide[read[pos_variant_read]&3];
                                        k++;
                                        pos_variant_read++;
                                }

                                newvar->alt[k] = '\0';
                                insert_variants(variant_list, newvar);
                        }
                } else if (code_align == CODE_DEL) {
                        int ps_var_genome = pos_variant_genome;
                        int ps_var_read = pos_variant_read;
                        variant_t *newvar = (variant_t*) malloc(sizeof(variant_t));

                        newvar->offset = ref_genome->pt_seq[A.coord_seq];
                        newvar->chr = ref_genome->seq_name[A.coord_seq];
                        newvar->depth = 1;
                        code_align_idx += 2;

                        while (code_align_tab[code_align_idx] < 4) {
                                ps_var_genome++;
                                code_align_idx++;
                        }

                        while (ref_genome->data[ps_var_genome] == read[ps_var_read]) {
                                ps_var_read--;
                                ps_var_genome--;
                                pos_variant_genome--;
                                pos_variant_read--;
                        }

                        newvar->pos = pos_variant_genome;
                        newvar->alt[0] = nucleotide[ref_genome->data[pos_variant_genome]];
                        newvar->alt[1] = '\0';

                        {
                                int k = 0;
                                while (pos_variant_genome <= ps_var_genome) {
                                        newvar->ref[k] = nucleotide[ref_genome->data[pos_variant_genome]&3];
                                        k++;
                                        pos_variant_genome++;
                                }

                                newvar->ref[k] = '\0';
                                insert_variants(variant_list, newvar);
                        }
                }
        }
}

int process_read(genome_t *ref_genome,
                 int8_t *reads_buffer,
                 variant_tree_t **variant_list,
                 int *substitution_list,
                 int8_t *mapping_coverage,
                 FILE *f1,
                 FILE *f2,
                 int round,
                 times_ctx_t *times_ctx,
                 reads_info_t *reads_info)
{
        double t1, t2;
        align_t *align_tab = (align_t *) malloc(sizeof(align_t) * MAX_ALIGN * NB_DPU);
        int nb_match = 0;
        char nucleotide[4] = {'A','C','T','G'};
        int nb_read_map = 0;
        int offset[4]; /* offset in align_tab of the first read  */
        int score[4];  /* min score                              */
        int nbread[4]; /* number of reads                        */
        int *align_tab_tmp = (int *) malloc(sizeof(int) * MAX_ALIGN);
        int size_read = reads_info->size_read;

        t1 = my_clock();

        /* Get output data from DPUs */
        for (int numdpu = 0; numdpu < NB_DPU; numdpu++) {
                int k = 0;
                int num_read = read_out_num(numdpu, k);

                while (num_read != -1) {
                        long coord = read_out_coord(numdpu, k);

                        align_tab[nb_match].num_read = num_read;
                        align_tab[nb_match].coord_seq = (int) (coord >> 32);
                        align_tab[nb_match].coord_pos = (int) (coord & 0xFFFFFFFF);
                        align_tab[nb_match].score = read_out_score(numdpu, k);
                        nb_match++;

                        num_read = read_out_num(numdpu, ++k);
                }
        }

        /* Sort output data from DPUs */
        qsort(align_tab, nb_match, sizeof(align_t), cmpalign);

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
                int align_idx = 0;
                int numpair = align_tab[i].num_read / 4;
                for (int k = 0; k < 4; k++) {
                        score[k] = 1000;
                        offset[k] = -1;
                        nbread[k] = 0;
                }
                while ((i < nb_match) && (numpair == align_tab[i].num_read / 4)) {
                        int type = align_tab[i].num_read % 4;
                        if (align_tab[i].score < score[type]) {
                                score[type] = align_tab[i].score;
                                offset[type] = i;
                                nbread[type] = 1;
                        } else {
                                if (align_tab[i].score == score[type]) nbread[type]++;
                        }
                        i++;
                }
                /* Compute a [0, 3] read pair */
                for (int pa = offset[0]; pa < offset[0] + nbread[0]; pa++) {
                        for (int pb = offset[3]; pb < offset[3] + nbread[3]; pb++) {
                                if (check_pair(align_tab[pa], align_tab[pb], reads_info)) {
                                        if (align_idx < MAX_ALIGN - 2) {
                                                align_tab_tmp[align_idx++] = pa;
                                                align_tab_tmp[align_idx++] = pb;
                                        }
                                }
                        }
                }
                /* Compute a [1, 2] read pair */
                for (int pa = offset[1]; pa < offset[1] + nbread[1]; pa++) {
                        for (int pb = offset[2]; pb < offset[2] + nbread[2]; pb++) {
                                if (check_pair(align_tab[pb], align_tab[pa], reads_info)) {
                                        if (align_idx < MAX_ALIGN - 2) {
                                                align_tab_tmp[align_idx++] = pa;
                                                align_tab_tmp[align_idx++] = pb;
                                        }
                                }
                        }
                }
                if (align_idx == 2) { /* Mapped reads*/
                        set_variant(align_tab[align_tab_tmp[0]],
                                    ref_genome,
                                    reads_buffer,
                                    variant_list,
                                    substitution_list,
                                    mapping_coverage,
                                    reads_info);
                        set_variant(align_tab[align_tab_tmp[1]],
                                    ref_genome,
                                    reads_buffer,
                                    variant_list,
                                    substitution_list,
                                    mapping_coverage,
                                    reads_info);
                        nb_read_map += 2;
                } else { /* Non-mapped reads.
                          * Add same to the reads to look at during next round
                          * (keep the already computed offset for next round).
                          */
                        int8_t *read = &reads_buffer[numpair * size_read];
                        fprintf(f1, ">>%d\n", SIZE_SEED * (round + 1));
                        for (int j = SIZE_SEED; j < size_read; j++) {
                                fprintf(f1, "%c", nucleotide[read[j]&3]);
                        }
                        for (int j = 0; j < SIZE_SEED; j++) {
                                fprintf(f1, "A");
                        }
                        fprintf(f1, "\n");
                        read = &reads_buffer[(numpair+2) * size_read];
                        fprintf(f2, ">>%d\n", SIZE_SEED * (round + 1));
                        for (int j = SIZE_SEED; j < size_read; j++) {
                                fprintf(f2, "%c", nucleotide[read[j]&3]);
                        }
                        for (int j = 0; j < SIZE_SEED; j++) {
                                fprintf(f2, "A");
                        }
                        fprintf(f2, "\n");
                }
        }

        free(align_tab_tmp);
        free(align_tab);

        t2 = my_clock();
        times_ctx->process_read = t2 - t1;
        times_ctx->tot_process_read += t2 - t1;

        return nb_read_map;
}
