/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vartree.h"
#include "genome.h"
#include "upvc.h"

#define MAX_BUFFER_SIZE (1024)

static int print_variant_tree(variant_tree_t *variant_tree, int8_t *mapping_coverage, FILE *vcf_file, FILE *variant_file)
{
        variant_t *variant;
        int nb_variant = 0;

        if (variant_tree == NULL) {
                return 0;
        }

        nb_variant = print_variant_tree(variant_tree->left, mapping_coverage, vcf_file, variant_file);

        variant = variant_tree->vars;
        if (variant->depth > 6) {
                fprintf(vcf_file,
                        "%s\t%llu\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d\n",
                        variant->chr,
                        (unsigned long long)((uint64_t)(variant->pos + 1) - variant->offset),
                        variant->ref,
                        variant->alt,
                        variant->depth,
                        mapping_coverage[variant->pos]);
                nb_variant++;
        }

        if (variant->depth > 2) {
                fprintf(variant_file,
                        "%s\t%llu\t%s\t%s\ttDEPTH=%d\tCOV=%d\n",
                        variant->chr,
                        (unsigned long long)((uint64_t)(variant->pos + 1) - variant->offset),
                        variant->ref,
                        variant->alt,
                        variant->depth,
                        mapping_coverage[variant->pos]);
        }

        nb_variant += print_variant_tree(variant_tree->right, mapping_coverage, vcf_file, variant_file);

        return nb_variant;
}

void create_vcf(char *chromosome_name,
                genome_t *ref_genome,
                variant_tree_t **variant_list,
                int *substitution_list,
                int8_t *mapping_coverage,
                times_ctx_t *times_ctx)
{
        double t1, t2;
        FILE *vcf_file, *variant_file;
        char nucleotide[4] = {'A','C','T','G'};
        int nb_variant;
        char filename[MAX_BUFFER_SIZE];

        printf("Create VCF\n");
        t1 = my_clock();

        sprintf(filename, "%svars_upvc.vcf", chromosome_name);
        vcf_file = fopen(filename,"w");

        sprintf(filename, "%svars_upvc.var", chromosome_name);
        variant_file = fopen(filename,"w");

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
        fprintf(vcf_file, "##reference=%s\n", ref_genome->fasta_file_name);

        /* print the column names (fields are tab-delimited in VCF) */
        fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");


        /* ####### END OF HEADER ####### */

        /* insert substitution variant in the variant_t tree */
        /* for each sequence in the genome */
        for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
                uint64_t start_position = ref_genome->pt_seq[seq_number];
                /* for each position in the sequence */
                for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
                        /* get the substitution at the `start_position+seq_position`th position in the genome */
                        uint64_t current_position = start_position + seq_position;
                        int substitution = substitution_list[current_position];
                        for (int i = 0; i < 4; i++) {
                                /* get number of substitutions for the base nt[i] at that position */
                                int nb_substitution = (substitution >> (i * 8)) & 0xFF;
                                if (nb_substitution > 7) {
                                        variant_t *newvar = (variant_t*)malloc(sizeof(variant_t));
                                        newvar->offset = start_position;
                                        newvar->chr = ref_genome->seq_name[seq_number];
                                        newvar->depth = nb_substitution;
                                        newvar->pos = current_position;
                                        newvar->ref[0] = nucleotide[ref_genome->data[current_position] & 3];
                                        newvar->ref[1] = '\0';
                                        newvar->alt[0] = nucleotide[i];
                                        newvar->alt[1] = '\0';
                                        insert_variants(variant_list, newvar);
                                }
                        }
                }
        }

        nb_variant = print_variant_tree(*variant_list, mapping_coverage, vcf_file, variant_file);

        fclose(vcf_file);

        t2 = my_clock();
        times_ctx->vcf = t2-t1;

        printf(" - number of variants: %d\n", nb_variant);
        printf(" - time: %lf sec.\n", times_ctx->vcf);
}
