#include <stdio.h>
#include <stdint.h>
#include "upvc.h"

// organisation de la memoire d'un DPU
typedef struct {
    // index
    int8_t *neighbor_idx;       // tableau qui contient le voisinage (code) des graines du genome de reference
    long *coordinate;           // tableau qui contient le numeros de sequence et l'indice dans a sequence de la graine

    // pool de reads
    int8_t *neighbor_read;      // tableau qui contient voisinage (code) du read
    int *offset;                // indique pour chaque read l'adresse du 1er voisinage
    int *count;                 // indique pour chaque read le nombre de voisinages
    int *num;                   // indique le numero des reads

    // resultat
    int *out_score;             // resultat : score
    int *out_num;               // resultat : no du read ou on a trouve un match
    long *out_coord;            // resultat : coordonnees sur le genome ou le read a matche
} MEM_DPU;

extern MEM_DPU *MDPU;

static void print_byte_to_sym(uint8_t byte, unsigned int offset, FILE *out) {
    uint8_t as_byte = (uint8_t) ((byte >> offset) & 0x3);
    char as_char;
    switch (as_byte) {
        case CODE_A:
            as_char = 'A';
            break;

        case CODE_C:
            as_char = 'C';
            break;

        case CODE_G:
            as_char = 'G';
            break;

        default:
            as_char = 'T';
    }
    fprintf(out, "%c", as_char);
}

void print_neighbor_idx(int d, int offs, int nr_nbr, FILE *out) {
    int each_nbr;

    for (each_nbr = 0; each_nbr < nr_nbr; each_nbr++) {
        fprintf(out, "\t");
        unsigned int i;
        for (i = 0; i < SIZE_NBR; i++) {
            uint8_t this_byte = (uint8_t) MDPU[d].neighbor_idx[(offs + each_nbr) * SIZE_NBR + i];
            print_byte_to_sym(this_byte, 0, out);
            print_byte_to_sym(this_byte, 2, out);
            print_byte_to_sym(this_byte, 4, out);
            print_byte_to_sym(this_byte, 6, out);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

void print_coordinates(int d, int offs, int l, FILE *out) {
    int i;
    for (i = 0; i < l; i++) {
        fprintf(out, " %lu", MDPU[d].coordinate[offs + i]);
    }
    fprintf(out, "\n");
}

void print_index_seeds(INDEX_SEED **SEED, FILE *out) {
    INDEX_SEED *C;
    int i;

    for (i = 0; i < NB_SEED; i++) {
        fprintf(out, "SEED=%u\n", i);
        C = SEED[i];
        while (C != NULL) {
            fprintf(out, "\tPU %u @%u[%u]\n", C->num_dpu, C->offset, C->nb_nbr);
            print_neighbor_idx(C->num_dpu, C->offset, C->nb_nbr, out);
            fprintf(out, "\t");
            print_coordinates(C->num_dpu, C->offset, C->nb_nbr, out);
            C = C->next;
        }
    }
}
