/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genome.h"
#include "parse_args.h"
#include "upvc.h"
#include "vartree.h"

#include "common.h"

static variant_t **variant_list[MAX_SEQ_GEN] = { NULL };
static pthread_mutex_t mutex;

void variant_tree_init()
{
    genome_t *genome = genome_get();
    pthread_mutex_init(&mutex, NULL);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        variant_list[each_seq] = (variant_t **)calloc(genome->len_seq[each_seq], sizeof(variant_t *));
    }
}

void variant_tree_insert(variant_t *var, uint32_t seq_nr, uint32_t offset_in_chr, int8_t *read)
{
    pthread_mutex_lock(&mutex);
    variant_t **entry = &variant_list[seq_nr][offset_in_chr];
    variant_t *vars = *entry;
    while (vars != NULL) {
        if (!strcmp(vars->ref, var->ref) && !strcmp(vars->alt, var->alt)) {
            vars->depth++;
            vars->score += var->score;
            free(var);
            var = vars;
            goto end;
        }
        vars = vars->next;
    }
    var->next = *entry;
    *entry = var;

end:
    var->reads = realloc(var->reads, var->depth * SIZE_READ);
    assert(var->reads != NULL);
    memcpy(&var->reads[SIZE_READ * (var->depth - 1)], read, SIZE_READ);

    pthread_mutex_unlock(&mutex);
}

void variant_tree_free()
{
    genome_t *genome = genome_get();
    pthread_mutex_destroy(&mutex);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        for (unsigned int i = 0; i < genome->len_seq[each_seq]; i++) {
            variant_t *tmp = variant_list[each_seq][i];
            while (tmp != NULL) {
                variant_t *to_free = tmp;
                tmp = tmp->next;
                free(to_free);
            }
        }
        free(variant_list[each_seq]);
    }
}

typedef struct {
    uint32_t percentage;
    uint32_t score;
} depth_filter_t;

depth_filter_t sub_filter[] = {
    [3] = { 16, 19 },
    [4] = { 17, 19 },
    [5] = { 18, 19 },
    [6] = { 19, 19 },
    [7] = { 20, 19 },
    [8] = { 20, 19 },
    [9] = { 20, 20 },
    [10] = { 20, 21 },
    [11] = { 20, 21 },
    [12] = { 20, 21 },
    [13] = { 20, 22 },
    [14] = { 20, 22 },
    [15] = { 20, 22 },
    [16] = { 20, 22 },
    [17] = { 20, 22 },
    [18] = { 20, 22 },
    [19] = { 20, 23 },
    [20] = { 20, 23 },
    [21] = { 20, 23 },
    [22] = { 20, 24 },
};

depth_filter_t indel_filter[] = {
    [2] = { 11, 15 },
    [3] = { 12, 16 },
    [4] = { 13, 20 },
    [5] = { 13, 20 },
    [6] = { 15, 22 },
    [7] = { 16, 23 },
    [8] = { 17, 23 },
    [9] = { 20, 23 },
    [10] = { 20, 23 },
    [11] = { 20, 23 },
    [12] = { 20, 25 },
    [13] = { 20, 25 },
    [14] = { 20, 25 },
    [15] = { 20, 25 },
    [16] = { 30, 25 },
    [17] = { 30, 30 },
    [18] = { 30, 40 },
};

static const char nucleotide[4] = { 'A', 'C', 'T', 'G' };

static bool print_variant_tree(variant_t *var, uint32_t seq_nr, uint64_t seq_pos, genome_t *ref_genome, FILE *vcf_file)
{
    char *chr = ref_genome->seq_name[seq_nr];
    uint32_t cov = ref_genome->mapping_coverage[seq_pos + ref_genome->pt_seq[seq_nr]];
    uint32_t depth = var->depth;
    uint32_t score = var->score / depth;
    uint32_t percentage = 100;
    if (cov != 0) {
        percentage = depth * 100 / cov;
    }

    if (strlen(var->ref) == 1 && strlen(var->alt) == 1) { /* SUBSTITUTION */
        if (depth < 3) {
            return false;
        } else if (depth > 22) {
            depth = 22;
        }
        if (!(score <= sub_filter[depth].score && percentage >= sub_filter[depth].percentage)) {
            return false;
        }
    } else { /* INSERTION OR DELETION */
        if (depth < 2) {
            return false;
        } else if (depth > 18) {
            depth = 18;
        }
        if (!(score <= indel_filter[depth].score && percentage >= indel_filter[depth].percentage)) {
            return false;
        }
    }

    fprintf(vcf_file, "%s\t%lu\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d;SCORE=%d", chr, seq_pos, var->ref, var->alt, var->depth, cov,
        score);

    for (unsigned int each_read = 0; each_read < depth; each_read++) {
        int8_t *read = &var->reads[each_read * SIZE_READ];
        fprintf(vcf_file, ";READ%u=", each_read);
        for (unsigned int i = 0; i < SIZE_READ; i++) {
            fprintf(vcf_file, "%c", nucleotide[read[i]]);
        }
    }
    fprintf(vcf_file, "\n");

    return true;
}

void create_vcf()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);

    FILE *vcf_file;
    int nb_variant = 0;
    char filename[1024];
    genome_t *ref_genome = genome_get();

    sprintf(filename, "%s_upvc.vcf", get_input_path());
    vcf_file = fopen(filename, "w");

    /* ####### START OF HEADER ####### */

    /* print vcf version (required) */
    fprintf(vcf_file, "##fileformat=VCFv4.3\n");

    /* print source of VCF file (this program) */
    fprintf(vcf_file, "##source=UPVC %s\n", VERSION);

    /* get the file date */
    char filedate[10];
    time_t mytime = time(NULL);
    strftime(filedate, 100, "%Y%d%m", localtime(&mytime));
    fprintf(vcf_file, "##fileDate=%s\n", filedate);

    /* print reference genome file name */
    fprintf(vcf_file, "##reference=%s\n", get_input_fasta());

    /* print the column names (fields are tab-delimited in VCF) */
    fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    /* ####### END OF HEADER ####### */

    /* for each sequence in the genome */
    for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
        /* for each position in the sequence */
        for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
            variant_t *var = variant_list[seq_number][seq_position];
            while (var != NULL) {
                nb_variant += print_variant_tree(var, seq_number, seq_position, ref_genome, vcf_file) ? 1 : 0;
                var = var->next;
            }
        }
    }

    fclose(vcf_file);

    printf("\tnumber of variants: %d\n", nb_variant);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}
