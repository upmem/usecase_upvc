#include <stdio.h>
#include <stdlib.h>
#include "compare.h"
#include "upvc.h"

#define PQD_INIT_VAL (99)

#define PATH_SUBSTITUTION (0)
#define PATH_INSERTION (1)
#define PATH_DELETION (2)

static int min(int a, int b)
{
        return a < b ? a : b;
}

void display_DPD_matrix(int **Matrix, reads_info_t *reads_info)
{
        int size_neighbour = reads_info->size_neighbour_in_32bits_words;

        for (int i = 0; i < size_neighbour / 2; i++) {
                for (int j = 0; j < size_neighbour / 2; j++) {
                        printf("%3d ", Matrix[i][j]);
                }
                printf("\n");
        }
        printf("---\n");
        for (int i = size_neighbour / 2; i < size_neighbour + 1; i++) {
                for (int j = size_neighbour / 2; j < size_neighbour + 1; j++) {
                        printf("%3d ", Matrix[i][j]);
                }
                printf("\n");
        }
}

static void DPD_compute(int s1,
                        int s2,
                        int *Dij,
                        int Dijm,
                        int Dimj,
                        int Dimjm,
                        int *Pij,
                        int Pijm,
                        int *Qij,
                        int Qimj,
                        int *path)
{
        int min_QP, d;

        *Pij = min(Dijm + COST_GAPO, Pijm + COST_GAPE);
        *Qij = min(Dimj + COST_GAPO, Qimj + COST_GAPE);

        if (*Pij < *Qij) {
                min_QP =  *Pij;
                *path = PATH_INSERTION;
        } else {
                min_QP = *Qij;
                *path = PATH_DELETION;
        }
        d = Dimjm;
        if ( (s1 & 3) != (s2 & 3) ) {
                d += COST_SUB;
        }
        if (d < min_QP) {
                *Dij = d;
                *path = PATH_SUBSTITUTION;
        } else {
                *Dij = min_QP;
        }

}

int DPD(int8_t *s1, int8_t *s2, backtrack_t *backtrack, reads_info_t *reads_info)
{
        int matrix_size = reads_info->size_neighbour_in_32bits_words + 1;
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
                        DPD_compute(s1[i - 1], s2[j - 1],
                                    &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1],
                                    &P[i][j], P[i][j - 1],
                                    &Q[i][j], Q[i - 1][j],
                                    &path[i][j]);
                }
                Q[i][i + diagonal] = PQD_INIT_VAL;
                D[i][i + diagonal] = PQD_INIT_VAL;
        }
        for (int i = diagonal; i < matrix_size - diagonal; i++) {
                P[i][i - diagonal] = PQD_INIT_VAL;
                D[i][i - diagonal] = PQD_INIT_VAL;
                for (int j = i - diagonal + 1; j < i + diagonal; j++) {
                        DPD_compute(s1[i - 1], s2[j - 1],
                                    &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1],
                                    &P[i][j], P[i][j - 1],
                                    &Q[i][j], Q[i - 1][j],
                                    &path[i][j]);
                }
                Q[i][i + diagonal]=PQD_INIT_VAL;
                D[i][i + diagonal]=PQD_INIT_VAL;
        }

        for (int i = matrix_size - diagonal; i < matrix_size; i++) {
                P[i][i - diagonal] = PQD_INIT_VAL;
                D[i][i - diagonal] = PQD_INIT_VAL;
                for (int j = i - diagonal + 1; j < matrix_size; j++) {
                        DPD_compute(s1[i - 1], s2[j - 1],
                                    &D[i][j], D[i][j - 1], D[i - 1][j], D[i - 1][j - 1],
                                    &P[i][j], P[i][j - 1],
                                    &Q[i][j], Q[i - 1][j],
                                    &path[i][j]);
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
                        } else if (path[i][j] == PATH_DELETION){
                                i--;
                        } else {
                                fprintf(stderr, "Error during deletion of indels\n");
                                exit(255);
                        }
                }

                /* i>1 && j>1 conditions erased the INDELS at the beginning */
                backtrack[0].type = CODE_END;
                while ( (i > 1) && (j > 1) ) {
                        if (path[i][j] == PATH_SUBSTITUTION) {
                                i--;
                                j--;
                                if (D[i][j] != D[i-1][j-1]) {
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
                                } else if (path[i][j] == PATH_DELETION){
                                        i--;
                                        backtrack[align_distance].type = CODE_DEL;
                                        backtrack[align_distance].ix = i;
                                        backtrack[align_distance].jx = j;
                                        align_distance++;
                                } else {
                                        fprintf(stderr, "Error during DPD compute\n");
                                        exit(255);
                                }
                        }
                }
        }

        return align_distance;
}
