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

#define MAX_BUFFER_SIZE (1024)

static void print_variant_tree(
    variant_t *var, uint32_t seq_nr, uint64_t seq_pos, genome_t *ref_genome, FILE *vcf_file, FILE *variant_file)
{
    char *chr = ref_genome->seq_name[seq_nr];
    fprintf(vcf_file, "%s\t%lu\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d;SCORE=%d\n", chr, seq_pos, var->ref, var->alt, var->depth,
        ref_genome->mapping_coverage[seq_pos + ref_genome->pt_seq[seq_nr]], var->score);
}

void create_vcf()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);

    FILE *vcf_file, *variant_file;
    int nb_variant = 0;
    char filename[MAX_BUFFER_SIZE];
    genome_t *ref_genome = genome_get();

    sprintf(filename, "%svars_upvc.vcf", get_input_path());
    vcf_file = fopen(filename, "w");

    sprintf(filename, "%svars_upvc.var", get_input_path());
    variant_file = fopen(filename, "w");

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
                nb_variant++;
                print_variant_tree(var, seq_number, seq_position, ref_genome, vcf_file, variant_file);
                var = var->next;
            }
        }
    }

    fclose(vcf_file);

    printf("\tnumber of variants: %d\n", nb_variant);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}
