/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "genome.h"
#include "parse_args.h"
#include "upvc.h"
#include "vartree.h"
#include "vcf.h"

#define MAX_BUFFER_SIZE (1024)

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

    fprintf(vcf_file, "%s\t%lu\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d;SCORE=%d\n", chr, seq_pos, var->ref, var->alt, var->depth, cov,
        score);

    return true;
}

void create_vcf()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);

    FILE *vcf_file;
    int nb_variant = 0;
    char filename[MAX_BUFFER_SIZE];
    genome_t *ref_genome = genome_get();

    sprintf(filename, "%svars_upvc.vcf", get_input_path());
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

    variant_t ***variant_tree = variant_tree_get();
    /* for each sequence in the genome */
    for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
        /* for each position in the sequence */
        for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
            variant_t *var = variant_tree[seq_number][seq_position];
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
