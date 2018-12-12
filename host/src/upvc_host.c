/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "parse_args.h"
#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "processread.h"
#include "vartree.h"
#include "vcf.h"
#include "genome.h"
#include "upvc_dpu.h"
#include "upvc.h"
#include "mram_dpu.h"
#include "backends_functions.h"
#include "simu_backend.h"
#include "dpu_backend.h"

static int map_var_call(char *filename_prefix,
                        int round,
                        devices_t *devices,
                        genome_t *ref_genome,
                        index_seed_t **index_seed,
                        variant_tree_t **variant_list,
                        int *substitution_list,
                        int8_t *mapping_coverage,
                        reads_info_t *reads_info,
                        times_ctx_t *times_ctx,
                        backends_functions_t *backends_functions)
{
        char filename[1024];
        FILE *fipe1, *fipe2;  /* pair-end input file descriptor */
        FILE *fope1, *fope2;  /* pair-end output file descriptor */
        int8_t *reads_buffer;
        int nb_read;
        int nb_read_total = 0;
        int nb_read_map = 0;
        int nb_pass = 0;
        unsigned int nb_dpu = get_nb_dpu();

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
                dispatch_t dispatch;

                nb_read_total += nb_read;

                printf("Round %d / Pass %d\n", round, nb_pass);
                printf(" - get %d reads (%d)\n", nb_read / 2, nb_read_total / 2);
                printf(" - time to get reads      : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->get_reads,
                       times_ctx->tot_get_reads);

                if (DEBUG_PASS != -1 && DEBUG_PASS != nb_pass) {
                        nb_pass++;
                        continue;
                }

                dispatch = dispatch_read(index_seed,
                                         reads_buffer,
                                         nb_read,
                                         nb_dpu,
                                         times_ctx,
                                         reads_info,
                                         backends_functions);
                printf(" - time to dispatch reads : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->dispatch_read,
                       times_ctx->tot_dispatch_read);

                backends_functions->run_dpu(dispatch,
                                            devices,
                                            (DEBUG_NB_RUN != -1) ? ((DEBUG_NB_RUN + DEBUG_FIRST_RUN) * get_nb_dpus_per_run()) : nb_dpu,
                                            times_ctx,
                                            reads_info);
                printf(" - time to map reads      : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->map_read,
                       times_ctx->tot_map_read);

                dispatch_free(dispatch, nb_dpu);

                if (DEBUG_PASS != -1) {
                        break;
                }

                nb_read_map += process_read(ref_genome,
                                            reads_buffer,
                                            variant_list,
                                            substitution_list,
                                            mapping_coverage,
                                            fope1,
                                            fope2,
                                            round,
                                            nb_dpu,
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

static void reload_and_verify_mram_images(reads_info_t *reads_info)
{
        index_seed_t **index_seed;
        FILE *seed_file;
        unsigned int nb_dpu = get_nb_dpu();
        malloc_dpu(reads_info, nb_dpu);
        index_seed = reload_mram_images_and_seeds(reads_info);
        seed_file = fopen(SEED_FILE_LOG, "w");
        print_index_seeds(index_seed, seed_file, reads_info);
        fclose(seed_file);
        printf("Please check %s to verify that the indexing is OK\n", SEED_FILE_LOG);

        free_index(index_seed);
        free_dpu(nb_dpu);
}

static void load_index_save_genome(reads_info_t *reads_info, times_ctx_t *times_ctx, backends_functions_t *backends_functions)
{
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);
        index_seed_t **index_seed = index_genome(ref_genome,
                                                 get_nb_dpu(),
                                                 times_ctx,
                                                 reads_info,
                                                 backends_functions);
        save_index_seeds(index_seed);

        free_genome(ref_genome);
        free_index(index_seed);
}

static void do_mapping(backends_functions_t *backends_functions, reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        index_seed_t **index_seed;
        char *input_prefix = get_input_path();
        unsigned int nb_dpu = get_nb_dpu();
        variant_tree_t *variant_list = NULL;
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);

        int8_t *mapping_coverage = (int8_t *) calloc(sizeof(int8_t), ref_genome->fasta_file_size);
        int *substitution_list = (int *) calloc(sizeof(int), ref_genome->fasta_file_size);

        index_seed = backends_functions->get_index_seed(nb_dpu, ref_genome, reads_info, times_ctx, backends_functions);

        if ((DEBUG_NB_RUN != -1 && DEBUG_FIRST_RUN == -1)
            || (DEBUG_NB_RUN == 0)) {
                ERROR_EXIT(42, "DEBUG MACRO has not been well configured!");
        }

        devices_t *devices = backends_functions->init_devices(get_nb_dpus_per_run(), get_dpu_binary());

        for (int round = 0; round < 3; round++) {
                if (DEBUG_ROUND != -1 && DEBUG_ROUND != round) {
                        continue;
                }
                map_var_call(input_prefix, round, devices, ref_genome, index_seed,
                             &variant_list, substitution_list, mapping_coverage,
                             reads_info, times_ctx, backends_functions);
        }

        backends_functions->free_devices(devices);

        create_vcf(input_prefix, ref_genome, &variant_list, substitution_list, mapping_coverage, times_ctx);

        free_variant_tree(variant_list);
        free_genome(ref_genome);
        free_index(index_seed);
        free(substitution_list);
        free(mapping_coverage);
        free_dpu(nb_dpu);
}

static void print_time()
{
        time_t timer;
        char time_buf[26];
        struct tm* tm_info;

        time(&timer);
        tm_info = localtime(&timer);

        strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
        printf("upvc started at: %s\n", time_buf);
}

int main(int argc, char *argv[])
{
        reads_info_t reads_info;
        times_ctx_t times_ctx;
        backends_functions_t backends_functions;

        memset(&times_ctx, 0, sizeof(times_ctx_t));

        validate_args(argc, argv);

        printf("%s\n", VERSION);
        print_time();

        reads_info.size_read = get_read_size(get_input_pe1());
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        printf("Information\n");
        printf(" - read size: %d\n", reads_info.size_read);

        setup_dpus_for_target_type(get_target_type());

        if (get_simulation_mode()) {
                backends_functions.init_devices = init_devices_simulation;
                backends_functions.free_devices = free_devices_simulation;
                backends_functions.run_dpu = run_dpu_simulation;
                backends_functions.add_seed_to_requests = add_seed_to_simulation_requests;
                backends_functions.get_index_seed = get_index_seed_simulation;
                backends_functions.init_vmis = init_vmis_simulation;
                backends_functions.free_vmis = free_vmis_simulation;
                backends_functions.write_vmi = write_vmi_simulation;
        } else {
                backends_functions.init_devices = dpu_try_alloc_for;
                backends_functions.free_devices = dpu_try_free;
                backends_functions.run_dpu = run_on_dpu;
                backends_functions.add_seed_to_requests = add_seed_to_dpu_requests;
                backends_functions.get_index_seed = get_index_seed_dpu;
                backends_functions.init_vmis = init_vmis_dpu;
                backends_functions.free_vmis = free_vmis_dpu;
                backends_functions.write_vmi = write_vmi_dpu;
        }

        switch(get_goal()) {
        case goal_index:
                load_index_save_genome(&reads_info, &times_ctx, &backends_functions);
                break;
        case goal_check:
                reload_and_verify_mram_images(&reads_info);
                break;
        case goal_map:
                do_mapping(&backends_functions, &reads_info, &times_ctx);
                break;
        case goal_unknown:
        default:
                ERROR_EXIT(23, "goal has not been specified!");
        }

        free_args();

        return 0;
}
