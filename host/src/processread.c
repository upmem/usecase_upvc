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
#include "mapping_file.h"
#include "upvc.h"
#include "vartree.h"
#include "parse_args.h"
#include "profiling.h"

#define DEBUG_READ_MAPPING false

#define SIZE_INSERT_MEAN (400)
#define SIZE_INSERT_STD (3 * 50)

#define CODE_A 0 /* ('A'>>1)&3   41H  0100 0001 */
#define CODE_C 1 /* ('C'>>1)&3   43H  0100 0011 */
#define CODE_T 2 /* ('T'>>1)&3   54H  0101 0100 */
#define CODE_G 3 /* ('G'>>1)&3   47H  0100 0111 */

#define PQD_INIT_VAL (999)

#if SIZE_READ>120
#define MAX_SUBSTITUTION 20
#else
#define MAX_SUBSTITUTION 20
#endif

#define MAX_SCORE_DIFFERENCE_WITH_BEST 40
#define MAX_CONSIDERED_MAPPINGS 4000


static int min(int a, int b) { return a < b ? a : b; }

static unsigned int reads_not_correct_cost = 0;
static unsigned int reads_correct_cost_sub_only = 0;
static unsigned int reads_correct_cost_DPD = 0;

static void DPD_compute(
    int s1, int s2, int *Dij, int Dijm, int Dimj, int Dimjm, int *Pij, int Pijm, int *Qij, int Qimj, int *xij)
{
    /* Compute values for D,P,Q and x at (i,j) from values at (i-1,j), (i,j-1), (i-1,j-1).
     * i represents relative index in reference genome (starting from where the read is mapped).
     * j represents index in read.
     * Dij is the minimum cost to align the first i nucleotides from genome with the first j nucleotides from read.
     * Pij is the same minimum cost but assuming the last operation was an insertion.
     * Qij is the same minimum cost but assuming the last operation was a deletion.
     * xij is the chosen operation to reach the cost in Dij :
     *      0 means read and reference genome match
     *      1 is for a substitution
     *      2 is for an insertion
     *      3 is for a deletion
     */
    int min_QP, d;

    // Compute cost in case of insertion :
    // D_i_j-1 + COST_GAPO in case of first insertion
    // D_i_j-1 + COST_GAPE in case insertion directly follows another insertion (In this case, D_i_j-1 = P_i_j-1)
    *Pij = min(Dijm + COST_GAPO, Pijm + COST_GAPE);
    // Similar for deletion :
    *Qij = min(Dimj + COST_GAPO, Qimj + COST_GAPE);
    *xij = 0;

    //x_i_j is the backtracking
    *xij = 0;

    int x;
    // Get minimum between Pij (which assumes insertion) and Qij (which assumes deletion)
    // and store the associated operation in x (2 for insertion and 3 for deletion).
    if (*Pij < *Qij) {
        min_QP = *Pij;
        x = 2;
    } else {
        min_QP = *Qij;
        x = 3;
    }

    // Compute the diagonal score in d :
    // Dimjm (D[i-1][j-1]) if genome and read match (genome[read_map_address+SIZE_SEED+i]==read[j]) 
    // Dimjm + COST_SUB if genome and read do not match (genome[read_map_address+SIZE_SEED+i]!=read[j])
    d = Dimjm;
    if ((s1 & 3) != (s2 & 3)) {
        d += COST_SUB;
        *xij = 1;
    }

    // If diagonal score is best, store it in Dij (in this case, xij has already been set to the correct value)
    // Otherwise, set Dij to the best score between insertion and deletion (and set xij correspondingly).
    if (d < min_QP) {
        *Dij = d;
    } else {
        *Dij = min_QP;
        *xij = x;
    }
}

int subOnlyPath(int8_t *s1, int8_t *s2, backtrack_t* backtrack, backtrack_t ** backtrack_end)
{
    STAT_RECORD_START(STAT_SUB_ONLY_PATH);
	int score = 0;
	backtrack[0].type = CODE_END;
	(*backtrack_end) = &backtrack[1];
	for (int i=0; i<SIZE_READ; i++) {
		if (s1[i] != s2[i]) {
                    (*backtrack_end)->type = CODE_SUB;
                    (*backtrack_end)->ix = i;
                    (*backtrack_end)->jx = i;
                    (*backtrack_end)++;
                    score += COST_SUB;
		} else {
                    (*backtrack_end)->type = 0;
                    (*backtrack_end)->ix = i;
                    (*backtrack_end)->jx = i;
                    (*backtrack_end)++;
		}
	}
    (*backtrack_end)--;
    STAT_RECORD_LAST_STEP(STAT_SUB_ONLY_PATH, 0);
	return score;
}

int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, backtrack_t ** backtrack_end)
{
    STAT_RECORD_START(STAT_DPD);
    /* s1 is a pointer to the genome where the end of the seed was mapped
     * s2 is a pointer to the end of the seed in the read
     * backtrack is a data structure that will be populated by this function with a sequence of operations
     *  (insertions, deletions, substitutions or match) paired with their corresponding indices in genome and read
     */
    // Matrices are only computed for neighbours (ie: nucleotides not in the seed).
    int matrix_size = SIZE_READ + 1;
    // Only NB_DIAG (odd) diagonals are computed for each matrix. Values outside this diagonal aren't considered.
    // ie: we only consider indices where j-diagonal < i < j+diagonal
    int diagonal = (NB_DIAG / 2) + 1;
    // D is the matrix of scores.
    // P is the matrix of scores assuming last operation is an insertion.
    // Q is the matrix of scores assuming last operation is a deletion.
    // X is the matrix of actual operations used for backtracking.
    int D[matrix_size][matrix_size];
    int P[matrix_size][matrix_size];
    int Q[matrix_size][matrix_size];
    int X[matrix_size][matrix_size];
    // Best score found at the end of computation (upper row/right-most column of D matrix)
    int min_score = PQD_INIT_VAL;
    // Indices of best score found in the upper-right row/column of D matrix
    int min_score_i_idx = 0;
    int min_score_j_idx = 0;


    // Set first row and column of P and Q to high values
    // Set first row and column of D to correct values
    // FIXME : It would probably make more sense to set D[i][0] and D[0][i] to i*COST_GAPO
    //         but this should also be changed accordingly in DPU then.
    for (int i = 0; i <= diagonal; i++) {
        P[i][0] = PQD_INIT_VAL;
        P[0][i] = PQD_INIT_VAL;
        Q[i][0] = PQD_INIT_VAL;
        Q[0][i] = PQD_INIT_VAL;
        D[i][0] = i * COST_SUB;
        D[0][i] = i * COST_SUB;
    }

    STAT_RECORD_STEP(STAT_DPD, 0);
    /* D matrix at this point : (assuming COST_SUB=1 and diagonal = 3 ie: NB_DIAG=5)
     * (0,0) in bottom left
     * 
     *  |
     *  |
     *  |3
     * ^|2
     * ||1
     * j|0 1 2 3
     *  +--------------
     *   i->
     */
    

    // Compute first trapezoid (bottom left part)
    for (int i = 1; i < diagonal; i++) {
        for (int j = 1; j < i + diagonal; j++) {
            DPD_compute(s1[i - 1], s2[j - 1], &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1], &P[i][j], P[i][j - 1],
                &Q[i][j], Q[i - 1][j], &X[i][j]);
        }
        Q[i][i + diagonal] = PQD_INIT_VAL;
        D[i][i + diagonal] = PQD_INIT_VAL;
    }
    STAT_RECORD_STEP(STAT_DPD, 1);

    /* D matrix at this point : (assuming COST_SUB=1 and diagonal = 3 ie: NB_DIAG=5)
     * (0,0) in bottom left
     * An x corresponds to a value that has been computed.
     * An M corresponds to PQD_INIT_VAL
     * 
     *  |
     *  |
     *  |    M
     *  |  M x
     *  |3 x x
     * ^|2 x x
     * ||1 x x
     * j|0 1 2 3
     *  +--------------
     *   i->
     */

    // Compute most of the diagonal
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
    STAT_RECORD_STEP(STAT_DPD, 2);

    /* D matrix at this point : (assuming COST_SUB=1, diagonal = 3 ie: NB_DIAG=5, and matrix_size=12)
     * (0,0) in bottom left
     * An x corresponds to a value that has been computed.
     * An M corresponds to PQD_INIT_VAL
     *  
     *  +-----------------------+
     *  |                M      |
     *  |              M x      |
     *  |            M x x      |
     *  |          M x x x      |
     *  |        M x x x x      |
     *  |      M x x x x x      |
     *  |    M x x x x x M      |
     *  |  M x x x x x M        |
     *  |3 x x x x x M          |
     * ^|2 x x x x M            |
     * ||1 x x x M              |
     * j|0 1 2 M                |
     *  +-----------------------+
     *   i->
     */

    // Compute last trapezoid (top right part)
    // (And check for best score in top-most row of D matrix)
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
    STAT_RECORD_STEP(STAT_DPD, 3);

    /* D matrix at this point : (assuming COST_SUB=1, diagonal = 3 ie: NB_DIAG=5, and matrix_size=12)
     * (0,0) in bottom left
     * An x corresponds to a value that has been computed.
     * An M corresponds to PQD_INIT_VAL
     *  
     *  +-----------------------+
     *  |                M x x x|
     *  |              M x x x x|
     *  |            M x x x x x|
     *  |          M x x x x x M|
     *  |        M x x x x x M  |
     *  |      M x x x x x M    |
     *  |    M x x x x x M      |
     *  |  M x x x x x M        |
     *  |3 x x x x x M          |
     * ^|2 x x x x M            |
     * ||1 x x x M              |
     * j|0 1 2 M                |
     *  +-----------------------+
     *   i->
     */

    // Find the best score in right-most column of D matrix (if it is better than the best from the top row)
    for (int j = matrix_size - diagonal; j < matrix_size; j++) {
        if (D[matrix_size - 1][j] < min_score) {
            min_score = D[matrix_size - 1][j];
            min_score_i_idx = matrix_size - 1;
            min_score_j_idx = j;
        }
    }
    STAT_RECORD_STEP(STAT_DPD, 4);

    // backtrack step
    {
        int i = min_score_i_idx;
        int j = min_score_j_idx;
        backtrack[0].type = CODE_END;
        (*backtrack_end) = &backtrack[1];
        while ((i > 0) || (j > 0)) {
            if(X[i][j] == 0) {
		        // Operation 0 : sequences match and nothing was done : decrease both indices.
		        i--;
		        j--;
		        (*backtrack_end)->type = CODE_MATCH;
		        (*backtrack_end)->ix = i;
		        (*backtrack_end)->jx = j;
		        (*backtrack_end)++;
	        } else if(X[i][j] == 1) {
		        // Operation 1 : substitution : decrease both indices and store substitution code.
		        i--;
		        j--;
		        (*backtrack_end)->type = CODE_SUB;
		        (*backtrack_end)->ix = i;
		        (*backtrack_end)->jx = j;
		        (*backtrack_end)++;
	        } else if(X[i][j] == 2) {
		        // Operation 2 : insertion : decrease read index but not genome index and store insertion code.
		        j--;
		        (*backtrack_end)->type = CODE_INS;
		        (*backtrack_end)->ix = i;
		        (*backtrack_end)->jx = j;
		        (*backtrack_end)++;
	        } else {
		        // Operation 3 : deletion : decrease genome index but not read index and store deletion code.
		        i--;
		        (*backtrack_end)->type = CODE_DEL;
		        (*backtrack_end)->ix = i;
		        (*backtrack_end)->jx = j;
		        (*backtrack_end)++;
	        }
	    }
        (*backtrack_end)--;
    }

    STAT_RECORD_LAST_STEP(STAT_DPD, 5);
    return min_score;
}

#define NB_FREQ_TABLE_MUTEXES (1<<10)
static pthread_mutex_t freq_table_mutexes[NB_FREQ_TABLE_MUTEXES+1];
//static pthread_mutex_t print_mutex;
static pthread_mutex_t codependence_list_mutex;

bool update_frequency_table(
		genome_t *ref_genome,
		dpu_result_out_t *result_tab,
		int8_t *reads_buffer,
		float *reads_quality_buffer,
		int pos,
		float mapq,
        bool unsure
		)
{
    STAT_RECORD_START(STAT_UPDATE_FREQUENCY_TABLE);
	struct frequency_info ** frequency_table = get_frequency_table();
    struct variants_codependence_info_list** codependence_table = get_codependence_table(&codependence_list_mutex);
	uint64_t genome_pos = ref_genome->pt_seq[result_tab[pos].coord.seq_nr] + result_tab[pos].coord.seed_nr;
	int num = result_tab[pos].num;
	int8_t *read = reads_buffer + (num * SIZE_READ);
	float *read_quality = reads_quality_buffer + (num>>1)*SIZE_READ;
	bool invert_read = num & 1;

	backtrack_t backtrack[SIZE_READ<<1];
	backtrack_t * backtrack_end = backtrack;

	/* First, try substitutions only */
	unsigned int computed_score = subOnlyPath(&ref_genome->data[genome_pos], read, backtrack, &backtrack_end);
    STAT_RECORD_STEP(STAT_UPDATE_FREQUENCY_TABLE, 0);

	if (computed_score > result_tab[pos].score) {
            computed_score = DPD(&ref_genome->data[genome_pos], read, backtrack, &backtrack_end);
            if (computed_score == result_tab[pos].score) {
                reads_correct_cost_DPD++;
            } else {
                reads_not_correct_cost++;
            }
	} else {
        reads_correct_cost_sub_only++;
    }
    STAT_RECORD_STEP(STAT_UPDATE_FREQUENCY_TABLE, 1);

    #if DEBUG_READ_MAPPING
    write_read_mapping_from_backtrack(ref_genome->seq_name[result_tab[pos].coord.seq_nr], result_tab[pos].coord.seed_nr, backtrack_end, read, num);
    #endif

    STAT_RECORD_STEP(STAT_UPDATE_FREQUENCY_TABLE, 2);

	// /!\ Pointer arithmetics
	unsigned int mutex_index = (genome_pos*NB_FREQ_TABLE_MUTEXES)/ref_genome->fasta_file_size;
	unsigned int end_mutex_index = ((genome_pos+backtrack_end->ix)*NB_FREQ_TABLE_MUTEXES)/ref_genome->fasta_file_size;
	pthread_mutex_lock(&freq_table_mutexes[mutex_index]);
	if (end_mutex_index != mutex_index) {
		pthread_mutex_lock(&freq_table_mutexes[end_mutex_index]);
	}
    STAT_RECORD_STEP(STAT_UPDATE_FREQUENCY_TABLE, 3);
	bool has_indel = false;
        uint8_t read_letter;
        //printf("genome_pos=%lu\n", genome_pos);
    unsigned int variants_found = 0;
    unsigned int variants_found_positions[SIZE_READ];
    uint8_t variants_found_letters[SIZE_READ];
	for (; backtrack_end > backtrack; backtrack_end--) {
            unsigned int current_position = genome_pos+backtrack_end->ix;
            //printf("backtrack_end->ix=%u; current_position=%u\n", backtrack_end->ix, current_position);
            if (current_position > ref_genome->fasta_file_size) {
                    continue;
            }
            switch (backtrack_end->type) {
                    case CODE_SUB:
                        read_letter = read[backtrack_end->jx];
                        if (read_letter > 3) {
                            read_letter = read_letter>>1 & 0x3;
                        }

                        // update codependence info
                        for (unsigned int i=0; i<variants_found; i++) {
                            //pthread_mutex_lock(&codependence_list_mutex);
                            add_codependence_info(&codependence_table[current_position], (int16_t) (variants_found_positions[i]-current_position), read_letter, variants_found_letters[i], ref_genome->fasta_file_size, &codependence_list_mutex);
                            add_codependence_info(&codependence_table[variants_found_positions[i]], (int16_t) (current_position-variants_found_positions[i]), variants_found_letters[i], read_letter, ref_genome->fasta_file_size, &codependence_list_mutex);
                            //pthread_mutex_unlock(&codependence_list_mutex);
                        }
                        variants_found_positions[variants_found] = current_position;
                        variants_found_letters[variants_found] = read_letter;
                        variants_found++;

                        frequency_table[read_letter][current_position].freq += mapq * read_quality[invert_read ? SIZE_READ-backtrack_end->jx-1 : backtrack_end->jx];
                        frequency_table[read_letter][current_position].score++;


                        if (unsure) {
                            frequency_table[read_letter][current_position].unsure_score++;
                        }
                        break;
                    case CODE_MATCH:
                        read_letter = read[backtrack_end->jx];
                        if (read_letter > 3) {
                            read_letter = read_letter>>1 & 0x3;
                        }
                        frequency_table[read_letter][current_position].freq += mapq * read_quality[invert_read ? SIZE_READ-backtrack_end->jx-1 : backtrack_end->jx];
                        frequency_table[read_letter][current_position].score++;


                        if (unsure) {
                            frequency_table[read_letter][current_position].unsure_score++;
                        }
                        break;
                    case CODE_INS:
                    case CODE_DEL:
                        has_indel = true;
                        //LOG_WARN("unhandled indel\n");
                        // TODO: handle indels
                        break;
            }
    }
    STAT_RECORD_STEP(STAT_UPDATE_FREQUENCY_TABLE, 4);
	pthread_mutex_unlock(&freq_table_mutexes[mutex_index]);
	if (end_mutex_index != mutex_index) {
		pthread_mutex_unlock(&freq_table_mutexes[end_mutex_index]);
	}
    STAT_RECORD_LAST_STEP(STAT_UPDATE_FREQUENCY_TABLE, 5);
	return has_indel;
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

static void set_variant(dpu_result_out_t result_match, genome_t *ref_genome, int8_t *reads_buffer)
{
    int8_t *read = &reads_buffer[result_match.num * SIZE_READ];
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    uint64_t genome_pos = ref_genome->pt_seq[result_match.coord.seq_nr] + result_match.coord.seed_nr;

    /* Get the differences betweend the read and the sequence of the reference genome that match */
    read = &reads_buffer[result_match.num * SIZE_READ];

    // code_alignment(code_result_tab, result_match.score, &ref_genome->data[genome_pos], read, size_neighbour_in_symbols);
	backtrack_t backtrack[SIZE_READ<<1];
	backtrack_t * backtrack_end = backtrack;

	/* First, try substitutions only */
	unsigned int computed_score = subOnlyPath(&ref_genome->data[genome_pos], read, backtrack, &backtrack_end);

	if (computed_score > result_match.score) {
            computed_score = DPD(&ref_genome->data[genome_pos], read, backtrack, &backtrack_end);
            if (computed_score == result_match.score) {
                reads_correct_cost_DPD++;
            } else {
                reads_not_correct_cost++;
            }
	} else {
        reads_correct_cost_sub_only++;
    }

    if (backtrack_end->type == CODE_ERR)
        return;

    /* Update "mapping_coverage" with the number of reads that match at this position of the genome */
    for (int i = 0; i < SIZE_READ; i++) {
        ref_genome->mapping_coverage[genome_pos + i] += 1;
    }

	for (; backtrack_end > backtrack; backtrack_end--) {
            unsigned int current_position = genome_pos+backtrack_end->ix;
            //printf("backtrack_end->ix=%u; current_position=%u\n", backtrack_end->ix, current_position);
            if (current_position > ref_genome->fasta_file_size) {
                    continue;
            }

            if (backtrack_end->type != 0) {
                variant_t *newvar = (variant_t *)malloc(sizeof(variant_t));
                newvar->depth = 1;
                newvar->score = result_match.score;
                newvar->next = NULL;
                int alt_idx = 0;
                int ref_idx = 0;
                int64_t pos_variant_genome = genome_pos + backtrack_end->ix;

                switch (backtrack_end->type) {
                    case CODE_INS:
                    case CODE_DEL:
                        pos_variant_genome--;
                        newvar->ref[ref_idx++] = nucleotide[ref_genome->data[genome_pos+(backtrack_end+1)->ix] & 0x3];
                        newvar->alt[alt_idx++] = nucleotide[read[(backtrack_end+1)->jx] & 0x3];
                        break;
                }

                for (;backtrack_end->type != 0 && backtrack_end->type != CODE_END; backtrack_end--) {
                    switch (backtrack_end->type) {
                        case CODE_SUB:
                            newvar->ref[ref_idx++] = nucleotide[ref_genome->data[genome_pos+backtrack_end->ix] & 0x3];
                            newvar->alt[alt_idx++] = nucleotide[read[backtrack_end->jx] & 0x3];
                            if (alt_idx >= MAX_SIZE_ALLELE || ref_idx >= MAX_SIZE_ALLELE) {
                                LOG_WARN("ignored read because of a too complex variant\n")
                            }
                            break;
                        case CODE_INS:
                            newvar->alt[alt_idx++] = nucleotide[read[backtrack_end->jx] & 0x3];
                            if (alt_idx >= MAX_SIZE_ALLELE) {
                                LOG_WARN("ignored read because of a too complex variant\n")
                            }
                            break;
                        case CODE_DEL:
                            newvar->ref[ref_idx++] = nucleotide[ref_genome->data[genome_pos+backtrack_end->ix] & 0x3];
                            if (ref_idx >= MAX_SIZE_ALLELE) {
                                LOG_WARN("ignored read because of a too complex variant\n")
                            }
                            break;
                    }
                }
                newvar->ref[ref_idx] = '\0';
                newvar->alt[alt_idx] = '\0';
                variant_tree_insert(
                                newvar, result_match.coord.seq_nr, pos_variant_genome - ref_genome->pt_seq[result_match.coord.seq_nr]);

            }
    }
}

static pthread_mutex_t non_mapped_mutex;
static void add_to_non_mapped_read(int numread, int round, FILE *fpe1, FILE *fpe2, int8_t *reads_buffer)
{
    STAT_RECORD_START(STAT_ADD_TO_NON_MAPPED_READ);
    if (fpe1 == NULL || fpe2 == NULL)
        return;
    pthread_mutex_lock(&non_mapped_mutex);
    STAT_RECORD_STEP(STAT_ADD_TO_NON_MAPPED_READ, 0);
    char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    int8_t *read = &reads_buffer[numread * SIZE_READ];
    fprintf(fpe1, ">>%d\n", SIZE_SEED * (round + 1));
    for (int j = SIZE_SEED; j < SIZE_READ; j++) {
        fprintf(fpe1, "%c", nucleotide[read[j] & 3]);
    }
    for (int j = 0; j < SIZE_SEED; j++) {
        fprintf(fpe1, "A");
    }
    fprintf(fpe1, "\n");
    read = &reads_buffer[(numread + 2) * SIZE_READ];
    fprintf(fpe2, ">>%d\n", SIZE_SEED * (round + 1));
    for (int j = SIZE_SEED; j < SIZE_READ; j++) {
        fprintf(fpe2, "%c", nucleotide[read[j] & 3]);
    }
    for (int j = 0; j < SIZE_SEED; j++) {
        fprintf(fpe2, "A");
    }
    fprintf(fpe2, "\n");
    STAT_RECORD_STEP(STAT_ADD_TO_NON_MAPPED_READ, 1);
    pthread_mutex_unlock(&non_mapped_mutex);
    STAT_RECORD_LAST_STEP(STAT_ADD_TO_NON_MAPPED_READ, 2);
}

/*#define USE_MAPQ_SCORE*/
static void do_process_read(process_read_arg_t *arg)
{
    STAT_RECORD_START(STAT_DO_PROCESS_READ);
    const unsigned int nb_match = arg->nb_match;
    dpu_result_out_t *result_tab = arg->result_tab;
    int round = arg->round;
    int8_t *reads_buffer = arg->reads_buffer;
    float *reads_quality_buffer = arg->reads_quality_buffer;
    genome_t *ref_genome = arg->ref_genome;
    FILE *fpe1 = arg->fpe1;
    FILE *fpe2 = arg->fpe2;

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
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 0);
      unsigned int i;
      if ((i = acquire_curr_match()) >= nb_match) {
        release_curr_match(i);
        STAT_RECORD_LAST_STEP(STAT_DO_PROCESS_READ, 6);
        return;
      }
      int numpair = result_tab[i].num / 4;
      unsigned int j = i;
      while ((j < nb_match) && (numpair == result_tab[j].num / 4)) {
          j++;
      }
      release_curr_match(j);
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 1);

      // find best mapping scores for each of the four reads considered
      unsigned int best_individual_scores[4] = {UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX};
      for (unsigned int x=i; x<j; x++)
      {
          unsigned int t = result_tab[x].num % 4;
          if (best_individual_scores[t] > result_tab[x].score)
              best_individual_scores[t] = result_tab[x].score;
      }

      unsigned int nb_considered_mappings[4] = {0,0,0,0};
      unsigned int considered_mappings[4][MAX_CONSIDERED_MAPPINGS];
      for (unsigned int x=i; x<j; x++)
      {
          int t = result_tab[x].num % 4;
          if (best_individual_scores[t] + MAX_SCORE_DIFFERENCE_WITH_BEST > result_tab[x].score && nb_considered_mappings[t]<MAX_CONSIDERED_MAPPINGS)
          {
              considered_mappings[t][nb_considered_mappings[t]++] = x;
          }
      }

      // i = start index in result_tab
      // j = stop index in result_tab
      // select best couples of paired reads
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 2);
      unsigned int P1[2] = {i,i};
      unsigned int P2[2] = {j-1,j-1};
      unsigned int pos1, pos2;
      unsigned int best_score[2] = { INVALID_SCORE, INVALID_SCORE };

      for (unsigned int x1 = 0; x1 < nb_considered_mappings[0]; x1++) {
        pos1 = result_tab[considered_mappings[0][x1]].coord.seed_nr;
        int score1 = result_tab[considered_mappings[0][x1]].score;
        for (unsigned int x2 = 0; x2 < nb_considered_mappings[3]; x2++) {
          pos2 = result_tab[considered_mappings[3][x2]].coord.seed_nr;
          int score2 = result_tab[considered_mappings[3][x2]].score;
          if ((abs((int)pos2 - (int)pos1) > READ_DIST_LOWER_BOUND && (abs((int)pos2 - (int)pos1) < READ_DIST_UPPER_BOUND))) {
            // update if this is one of the two best scores
            keep_best_2_scores(score1 + score2, P1, P2, considered_mappings[0][x1], considered_mappings[3][x2], best_score);
          }
        }
      }
      for (unsigned int x1 = 0; x1 < nb_considered_mappings[1]; x1++) {
        pos1 = result_tab[considered_mappings[1][x1]].coord.seed_nr;
        int score1 = result_tab[considered_mappings[1][x1]].score;
        for (unsigned int x2 = 0; x2 < nb_considered_mappings[2]; x2++) {
          pos2 = result_tab[considered_mappings[2][x2]].coord.seed_nr;
          int score2 = result_tab[considered_mappings[2][x2]].score;
          if ((abs((int)pos2 - (int)pos1) > READ_DIST_LOWER_BOUND && (abs((int)pos2 - (int)pos1) < READ_DIST_UPPER_BOUND))) {
            // update if this is one of the two best scores
            keep_best_2_scores(score1 + score2, P1, P2, considered_mappings[1][x1], considered_mappings[2][x2], best_score);
          }
        }
      }
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 3);

      bool update = false;
      bool hasIndel = false;

      unsigned np = get_nb_scores(best_score);
      if (np > 0) {
	LOG_DEBUG("found at least a pair (%u)\n", np);

        if(np == 2) {//if(false) {

          // found at least 2 matching pairs of positions. Check the delta between the two pairs to 
          // decide whether we should keep the best pair
          int delta = abs((int)(MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]])) 
              - (int)(MISMATCH_COUNT(result_tab[P1[1]]) + MISMATCH_COUNT(result_tab[P2[1]])));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) 
            + MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((2 * (MAX_SUBSTITUTION + 1)) - delta);
          if(delta_corrected < 0) {
            LOG_WARN("negative delta for square root %d\n", delta_corrected);
          }
          else if(delta >= DIST_PAIR_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta >= DIST_PAIR_THRESHOLD) {
#endif
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, false);
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, false);
            } else {
                set_variant(result_tab[P1[0]], ref_genome, reads_buffer);
                set_variant(result_tab[P2[0]], ref_genome, reads_buffer);
            }
            update = true;
          } else {
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, true);
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, true);
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[1], mapq, true);
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[1], mapq, true);
            } else {
                set_variant(result_tab[P1[0]], ref_genome, reads_buffer);
                set_variant(result_tab[P2[0]], ref_genome, reads_buffer);
                set_variant(result_tab[P1[1]], ref_genome, reads_buffer);
                set_variant(result_tab[P2[1]], ref_genome, reads_buffer);
            }
            // LOG_WARN("unusable pair (%u)\n", result_tab[i].num/4);
          }
        }
        else if(np) { // only one result, take it
          int delta = abs((int)(MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]])) - (2 * (MAX_SUBSTITUTION + 1)));
          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((2 * (MAX_SUBSTITUTION + 1)) - delta);
          if(delta_corrected < 0) {
            LOG_WARN("negative delta (np == 1) for square root %d\n", delta_corrected);
          }
          else if(delta >= DIST_PAIR_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta >= DIST_PAIR_THRESHOLD) {
#endif
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, false);
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, false);
            } else {
                set_variant(result_tab[P1[0]], ref_genome, reads_buffer);
                set_variant(result_tab[P2[0]], ref_genome, reads_buffer);
            }
            update = true;
          }/* else {
            LOG_WARN("unusable pair (%u)\n", result_tab[i].num/4);
          }*/
        }
      }// else {
      if (true) {

        // check mapping of R1 and R2 independently
        unsigned int best_score_R1[2] = { INVALID_SCORE, INVALID_SCORE };
        unsigned int best_score_R2[2] = { INVALID_SCORE, INVALID_SCORE };
        unsigned int throwaway[2];
        P1[0] = 0;
        P2[0] = 0;
        P1[1] = 0;
        P2[1] = 0;
        for (unsigned int read = i; read < j; read++) {
          unsigned t1 = result_tab[read].num % 4;
          if(t1 < 2) { // PE1 or RPE1
            keep_best_2_scores(result_tab[read].score, P1, throwaway, read, 0, best_score_R1);
          }
          else { // PE2 or RPE2
            keep_best_2_scores(result_tab[read].score, throwaway, P2, 0, read, best_score_R2);
          }
        }

        unsigned np1 = get_nb_scores(best_score_R1), np2 = get_nb_scores(best_score_R2);
        if(np1 == 2) {

          int delta = abs((int)MISMATCH_COUNT(result_tab[P1[0]]) - (int)MISMATCH_COUNT(result_tab[P1[1]]));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            LOG_WARN("negative delta (np1 == 2) for square root %d\n", delta_corrected);
          }
          else if(delta >= DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta >= DIST_SINGLE_THRESHOLD) {
#endif
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, false);
            } else {
                set_variant(result_tab[P1[0]], ref_genome, reads_buffer);
            }
            update = true;
          }
        }
        else if(np1) {
          int delta = abs((int)MISMATCH_COUNT(result_tab[P1[0]]) - (MAX_SUBSTITUTION + 1));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P1[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            LOG_WARN("negative delta (np1 == 1) for square root %d\n", delta_corrected);
          }
          else if(delta >= DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta >= DIST_SINGLE_THRESHOLD) {
#endif
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P1[0], mapq, false);
            } else {
                set_variant(result_tab[P1[0]], ref_genome, reads_buffer);
            }
            update = true;
          }
        }

        if(np2 == 2) {

          int delta = abs((int)MISMATCH_COUNT(result_tab[P2[0]]) - (int)MISMATCH_COUNT(result_tab[P2[1]]));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            LOG_WARN("negative delta (np2 == 2) for square root %d, %d %d %d %d\n", 
                delta_corrected, MISMATCH_COUNT(result_tab[P2[0]]), MISMATCH_COUNT(result_tab[P2[1]]), MAX_SUBSTITUTION + 1, delta);
          }
          else if(delta >= DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if(delta >= DIST_SINGLE_THRESHOLD) {
#endif

            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, false);
            } else {
                set_variant(result_tab[P2[0]], ref_genome, reads_buffer);
            }
            update = true;
          }
        }
        else if(np2) {
          int delta = abs((int)MISMATCH_COUNT(result_tab[P2[0]]) - (MAX_SUBSTITUTION + 1));

          float mapq = 1.0f;
#ifdef USE_MAPQ_SCORE
          int delta_corrected = MISMATCH_COUNT(result_tab[P2[0]]) + MAPQ_SCALING_FACTOR * ((MAX_SUBSTITUTION + 1) - delta);

          if(delta_corrected < 0) {
            LOG_WARN("negative delta (np2 == 1) for square root %d\n", delta_corrected);
          }
          else if (delta >= DIST_SINGLE_THRESHOLD) {
            mapq = 1.0 - sqrt((double)delta_corrected / SIZE_READ);
#else
          if (delta >= DIST_SINGLE_THRESHOLD) {
#endif
            if (get_use_frequency_table()) {
                hasIndel |= update_frequency_table(ref_genome, result_tab, reads_buffer, reads_quality_buffer, P2[0], mapq, false);
            } else {
                set_variant(result_tab[P2[0]], ref_genome, reads_buffer);
            }
            update = true;
          }
        }
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 4);
        if(!update) {
            //pthread_mutex_lock(&nr_reads_mutex);
            nr_reads_non_mapped++;
            //pthread_mutex_unlock(&nr_reads_mutex);
            add_to_non_mapped_read(numpair * 4, round, fpe1, fpe2, reads_buffer);
        }
        if(hasIndel)
            nr_reads_with_indels++;
        pthread_mutex_lock(&nr_reads_mutex);
        nr_reads_total_from_dpus++;
        pthread_mutex_unlock(&nr_reads_mutex);
      }
      STAT_RECORD_STEP(STAT_DO_PROCESS_READ, 5);
    }
    STAT_RECORD_LAST_STEP(STAT_DO_PROCESS_READ, 6);
}

#define PROCESS_READ_THREAD (8)
#define PROCESS_READ_THREAD_SLAVE (PROCESS_READ_THREAD - 1)
static process_read_arg_t args;
static pthread_t thread_id[PROCESS_READ_THREAD_SLAVE];
static bool stop_threads = false;

void process_read(FILE *fpe1, FILE *fpe2, int round, unsigned int pass_id)
{
    STAT_RECORD_START(STAT_PROCESS_READ);
    int8_t *reads_buffer = get_reads_buffer(pass_id);
    float *reads_quality_buffer = get_reads_quality_buffer(pass_id);
    acc_results_t acc_res = accumulate_get_result(pass_id, true);
    nr_reads_total += get_reads_in_buffer(pass_id) / 4;
    curr_match = 0;

    args.nb_match = acc_res.nb_res;
    args.result_tab = acc_res.results;
    args.round = round;
    args.reads_buffer = reads_buffer;
    args.reads_quality_buffer = reads_quality_buffer;
    args.fpe1 = fpe1;
    args.fpe2 = fpe2;

    STAT_RECORD_STEP(STAT_PROCESS_READ, 0);
    pthread_barrier_wait(&barrier);
    STAT_RECORD_STEP(STAT_PROCESS_READ, 1);
    do_process_read(&args);
    STAT_RECORD_STEP(STAT_PROCESS_READ, 2);
    pthread_barrier_wait(&barrier);
    STAT_RECORD_STEP(STAT_PROCESS_READ, 3);

    free(acc_res.results);
    STAT_RECORD_LAST_STEP(STAT_PROCESS_READ, 4);
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
    open_mapping_file();
#endif
    genome_t *ref_genome = genome_get();
    args.ref_genome = ref_genome;

    assert(pthread_mutex_init(&curr_match_mutex, NULL) == 0);
    assert(pthread_mutex_init(&non_mapped_mutex, NULL) == 0);
    for (int i=0; i<NB_FREQ_TABLE_MUTEXES+1; i++) {
        assert(pthread_mutex_init(&freq_table_mutexes[i], NULL) == 0);
    }
    assert(pthread_barrier_init(&barrier, NULL, PROCESS_READ_THREAD) == 0);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_create(&thread_id[each_thread], NULL, process_read_thread_fct, &args) == 0);
    }
}

void process_read_free()
{
#if DEBUG_READ_MAPPING
    close_mapping_file();
#endif
    stop_threads = true;
    pthread_barrier_wait(&barrier);

    for (unsigned int each_thread = 0; each_thread < PROCESS_READ_THREAD_SLAVE; each_thread++) {
        assert(pthread_join(thread_id[each_thread], NULL) == 0);
    }

    assert(pthread_barrier_destroy(&barrier) == 0);
    assert(pthread_mutex_destroy(&curr_match_mutex) == 0);
    assert(pthread_mutex_destroy(&non_mapped_mutex) == 0);
    for (int i=0; i<NB_FREQ_TABLE_MUTEXES+1; i++) {
        assert(pthread_mutex_init(&freq_table_mutexes[i], NULL) == 0);
    }
    assert(pthread_mutex_init(&codependence_list_mutex, NULL) == 0);
    printf("reads for which score not correctly computed: %u\n", reads_not_correct_cost);
    printf("reads for which score was correctly computed with DPD: %u\n", reads_correct_cost_DPD);
    printf("reads for which score was correctly computed with sub only: %u\n", reads_correct_cost_sub_only);
    fflush(stdout);
    fprintf(stderr, "%% reads non mapped: %f%%\n", (float)nr_reads_non_mapped * 100.0 / (float)nr_reads_total_from_dpus);
    fprintf(stderr, "%% reads with indels: %f%%\n", (float)nr_reads_with_indels * 100.0 / (float)(nr_reads_total_from_dpus - nr_reads_non_mapped));
    fprintf(stderr, "%% Total reads from dpus: %ld\n", nr_reads_total_from_dpus);
    fprintf(stderr, "%% Total reads: %ld\n", nr_reads_total);
}
