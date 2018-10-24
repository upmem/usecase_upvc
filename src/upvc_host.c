#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "upvc_dpu.h"
#include "processread.h"
#include "vartree.h"
#include "vcf.h"
#include "genome.h"
#include "upvc.h"

static void run_dpu(char *prog_to_offload, times_ctx_t *times_ctx, reads_info_t *reads_info)
{
        int numdpu;
        double t1, t2;
        pthread_t thread_id[NB_DPU];
        align_on_dpu_arg_t thread_args[NB_DPU];

        t1 = my_clock();

        for (numdpu = 0; numdpu < NB_DPU; numdpu++) {
                thread_args[numdpu].reads_info = *reads_info;
                thread_args[numdpu].numdpu = numdpu;
        }

        if (strcmp(prog_to_offload, "mapping") == 0) {
                for (numdpu = 0; numdpu < NB_DPU; numdpu++) {
                        pthread_create(&thread_id[numdpu], NULL, align_on_dpu, &thread_args[numdpu]);
                }
                for (numdpu = 0; numdpu < NB_DPU; numdpu++) {
                        pthread_join(thread_id[numdpu], NULL);
                }
        }

        t2 = my_clock();
        times_ctx->map_read = t2 - t1;
        times_ctx->tot_map_read += t2 - t1;
}

static int map_var_call(char *filename_prefix,
                        int round,
                        genome_t *ref_genome,
                        index_seed_t **index_seed,
                        variant_tree_t **variant_list,
                        int *substitution_list,
                        int8_t *mapping_coverage,
                        reads_info_t *reads_info,
                        times_ctx_t *times_ctx)
{
        char filename[1024];
        FILE *fipe1, *fipe2;  // pair-end input file descriptor
        FILE *fope1, *fope2;  // pair-end output file descriptor
        int8_t *reads_buffer;
        int nb_read;
        int nb_read_total = 0;
        int nb_read_map = 0;
        int nb_pass = 0;

        reads_info->delta_neighbour_in_bytes = (SIZE_SEED * round)/4;
        reads_info->size_neighbour_in_32bits_words =
                (reads_info->size_neighbour_in_bytes-reads_info->delta_neighbour_in_bytes) * 4;

        if (round == 0) {
                sprintf(filename, "%s_PE1.fastq", filename_prefix);
                fipe1 = fopen(filename, "r");
                sprintf(filename, "%s_PE2.fastq", filename_prefix);
                fipe2 = fopen(filename, "r");
        } else {
                sprintf(filename, "%s_%d_PE1.fasta", filename_prefix, round);
                fipe1 = fopen(filename, "r");
                sprintf(filename, "%s_%d_PE2.fasta", filename_prefix, round);
                fipe2 = fopen(filename, "r");
        }

        sprintf(filename, "%s_%d_PE1.fasta", filename_prefix, round+1);
        fope1 = fopen(filename, "w");
        sprintf(filename, "%s_%d_PE2.fasta", filename_prefix, round+1);
        fope2 = fopen(filename, "w");

        reads_buffer = (int8_t *) malloc(sizeof(int8_t) * MAX_READS_BUFFER * reads_info->size_read);

        /*
         * Loop:
         *   - Read a group of reads
         *   - Dispatch reads on DPUs
         *   - Execution on DPUs
         *   - Reads post-processing
         */
        while ((nb_read = get_reads(fipe1, fipe2, reads_buffer, times_ctx, reads_info)) != 0) {
                nb_read_total += nb_read;

                printf("Round %d / Pass %d\n", round, nb_pass);
                printf(" - get %d reads (%d)\n", nb_read / 2, nb_read_total / 2);
                printf(" - time to get reads      : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->get_reads,
                       times_ctx->tot_get_reads);

                dispatch_read(index_seed, reads_buffer, nb_read, times_ctx, reads_info);
                printf(" - time to dispatch reads : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->dispatch_read,
                       times_ctx->tot_dispatch_read);

                run_dpu("mapping", times_ctx, reads_info);
                printf(" - time to map reads      : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->map_read,
                       times_ctx->tot_map_read);

                nb_read_map += process_read(ref_genome,
                                            reads_buffer,
                                            variant_list,
                                            substitution_list,
                                            mapping_coverage,
                                            fope1,
                                            fope2,
                                            round,
                                            times_ctx,
                                            reads_info);
                printf(" - time to process reads  : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->process_read,
                       times_ctx->tot_process_read);
                printf(" - map %d reads\n", nb_read_map);
                printf("\n");
                nb_pass++;
        }

        free(reads_buffer);
        fclose(fipe1);
        fclose(fipe2);
        fclose(fope1);
        fclose(fope2);

        return nb_read_map;
}


int main(int argc, char *argv[])
{
        char filename[1024];
        genome_t *ref_genome;
        index_seed_t **index_seed;
        variant_tree_t *variant_list = NULL;
        int *substitution_list;
        int8_t *mapping_coverage;
        reads_info_t reads_info;
        times_ctx_t times_ctx = {0.0};

        if (argc != 2){
                printf("\nusage:\n  %s <genome> \n\n", argv[0]);
                exit(255);
        }
        printf("%s\n", VERSION);

        sprintf(filename, "%s.fasta", argv[1]);
        int accessref = access(filename, R_OK) == 0;
        sprintf(filename, "%s_PE1.fastq", argv[1]);
        int accesspe1 = access(filename, R_OK) == 0;
        sprintf(filename, "%s_PE2.fastq", argv[1]);
        int accesspe2 = access(filename, R_OK) == 0;

        if (accessref + accesspe1 + accesspe2 != 3) {
                printf("\nCould not find one or more of these files in the current directory:"
                       " %s.fasta, %s_PE1.fastq, %s_PE2.fastq \n\n",
                       argv[1], argv[1], argv[1]);
                exit(255);
        }

        reads_info.size_read = get_read_size(argv[1]);
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        printf("Information\n");
        printf(" - read size: %d\n", reads_info.size_read);

        printf("Read genome\n");
        ref_genome = get_genome(argv[1], &times_ctx);
        printf(" - #seq: %d\n", ref_genome->nb_seq);
        printf(" - time: %lf\n", times_ctx.get_genome);

        printf("Index genome\n");
        index_seed = index_genome(ref_genome, &times_ctx, &reads_info);
        printf(" - time: %lf\n", times_ctx.index_genome);

        mapping_coverage = (int8_t *) calloc(sizeof(int8_t), ref_genome->fasta_file_size);
        substitution_list = (int *) calloc(sizeof(int), ref_genome->fasta_file_size);


        map_var_call(argv[1], 0, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     &reads_info, &times_ctx);
        map_var_call(argv[1], 1, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     &reads_info, &times_ctx);
        map_var_call(argv[1], 2, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     &reads_info, &times_ctx);


        printf("Create VCF\n");
        int nb_var = create_vcf(argv[1], ref_genome, &variant_list, substitution_list, mapping_coverage, &times_ctx);
        printf(" - number of variants: %d\n", nb_var);
        printf(" - time: %lf sec.\n", times_ctx.vcf);

        free_variant_tree(variant_list);
        free_genome(ref_genome);
        free_index(index_seed);
        free(substitution_list);
        free(mapping_coverage);
        free_dpu();

        return 0;
}
