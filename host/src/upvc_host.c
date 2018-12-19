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

#define MAX_NB_PASS (1024)

static void run_pass(int round,
                     int nb_read,
                     unsigned int dpu_offset,
                     int nb_pass,
                     int8_t *reads_buffer,
                     dpu_result_out_t **result_tab,
                     unsigned int *result_tab_nb_read,
                     index_seed_t **index_seed,
                     dispatch_request_t *dispatch_requests,
                     devices_t *devices,
                     reads_info_t *reads_info,
                     times_ctx_t *times_ctx,
                     backends_functions_t *backends_functions)
{
        unsigned nb_read_map;
        printf("Round %d / DPU offset %d / Pass %d\n", round, dpu_offset, nb_pass);

        if (DEBUG_PASS != -1 && DEBUG_PASS != nb_pass) {
                return;
        }

        PRINT_TIME_DISPATCH(times_ctx, nb_pass);
        dispatch_read(index_seed,
                      reads_buffer,
                      nb_read,
                      dispatch_requests,
                      times_ctx,
                      reads_info,
                      backends_functions);
        printf(" - time to dispatch reads : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->dispatch_read,
               times_ctx->tot_dispatch_read);
        PRINT_TIME_DISPATCH(times_ctx, nb_pass);

        PRINT_TIME_MAP_READ(times_ctx, nb_pass);
        backends_functions->run_dpu(dispatch_requests,
                                    devices,
                                    dpu_offset,
                                    nb_pass,
                                    times_ctx,
                                    reads_info);
        printf(" - time to write reads      : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->write_reads,
               times_ctx->tot_write_reads);
        printf(" - time to compute          : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->compute,
               times_ctx->tot_compute);
        printf(" - time to read results     : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->read_result,
               times_ctx->tot_read_result);
        printf(" - time to map reads        : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->map_read,
               times_ctx->tot_map_read);
        PRINT_TIME_MAP_READ(times_ctx, nb_pass);


        if (DEBUG_PASS != -1) {
                return;
        }

        PRINT_TIME_ACC_READ(times_ctx, nb_pass);
        nb_read_map = accumulate_read(result_tab, result_tab_nb_read, dpu_offset, times_ctx);
        printf(" - time to accumulate read  : %7.2lf sec. / %7.2lf sec.\n",
               times_ctx->acc_read,
               times_ctx->tot_acc_read);
        printf(" - map %i reads\n", nb_read_map);
        printf("\n");
        PRINT_TIME_ACC_READ(times_ctx, nb_pass);
}

static void map_var_call(int round,
                         unsigned int nb_pass_total,
                         devices_t *devices,
                         genome_t *ref_genome,
                         index_seed_t **index_seed,
                         dispatch_request_t *dispatch_requests,
                         variant_tree_t **variant_list,
                         int *substitution_list,
                         int8_t *mapping_coverage,
                         int8_t **reads_buffer,
                         int *nb_read,
                         FILE *fope1,
                         FILE *fope2,
                         reads_info_t *reads_info,
                         times_ctx_t *times_ctx,
                         backends_functions_t *backends_functions)
{
        unsigned int nb_dpu = get_nb_dpu();
        unsigned int nb_read_map;
        unsigned int result_tab_nb_read[nb_pass_total];
        dpu_result_out_t *result_tab[nb_pass_total];

        memset(result_tab, 0, sizeof(dpu_result_out_t *) * nb_pass_total);
        memset(result_tab_nb_read, 0, sizeof(unsigned int) * nb_pass_total);

        /*
         * Loop:
         *   - Read a group of reads
         *   - Dispatch reads on DPUs
         *   - Execution on DPUs
         *   - Reads post-processing
         */
        for (unsigned int dpu_offset = 0; dpu_offset < nb_dpu; dpu_offset += get_nb_dpus_per_run()) {

                PRINT_TIME_WRITE_MRAM(times_ctx, 0);
                backends_functions->load_mram(dpu_offset, devices, reads_info, times_ctx);
                printf(" - time to write MRAMs : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->write_mram,
                       times_ctx->tot_write_mram);

                PRINT_TIME_WRITE_MRAM(times_ctx, 0);

                for (unsigned int each_pass = 0; each_pass < nb_pass_total; each_pass++) {
                        run_pass(round, nb_read[each_pass], dpu_offset, each_pass,
                                 reads_buffer[each_pass], &result_tab[each_pass], &result_tab_nb_read[each_pass],
                                 index_seed, dispatch_requests, devices,
                                 reads_info, times_ctx, backends_functions);
                }

        }

        for (unsigned int each_pass = 0; each_pass < nb_pass_total; each_pass++) {
                nb_read_map = process_read(ref_genome,
                                           reads_buffer[each_pass],
                                           variant_list,
                                           substitution_list,
                                           mapping_coverage,
                                           result_tab[each_pass],
                                           result_tab_nb_read[each_pass],
                                           fope1,
                                           fope2,
                                           round,
                                           times_ctx,
                                           reads_info);
                printf(" - time to process reads (pass %i) : %7.2lf sec. / %7.2lf sec.\n",
                       each_pass,
                       times_ctx->process_read,
                       times_ctx->tot_process_read);
                printf(" - map %i reads (pass %i)\n", nb_read_map, each_pass);
                free(reads_buffer[each_pass]);
        }
}

static void reload_and_verify_mram_images(reads_info_t *reads_info)
{
        index_seed_t **index_seed;
        FILE *seed_file;
        unsigned int nb_dpu = get_nb_dpu();
        malloc_dpu(reads_info, nb_dpu);
        index_seed = load_index_seeds();
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
        free_dpu(get_nb_dpu());
}

static void do_mapping(backends_functions_t *backends_functions, reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        index_seed_t **index_seed;
        char *input_prefix = get_input_path();
        unsigned int nb_dpu = get_nb_dpu();
        variant_tree_t *variant_list = NULL;
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);
        devices_t *devices;

        int8_t *mapping_coverage = (int8_t *) calloc(sizeof(int8_t), ref_genome->fasta_file_size);
        int *substitution_list = (int *) calloc(sizeof(int), ref_genome->fasta_file_size);
        dispatch_request_t *dispatch_requests = dispatch_create(nb_dpu, reads_info);

        if ((DEBUG_NB_RUN != -1 && DEBUG_FIRST_RUN == -1)
            || (DEBUG_NB_RUN == 0)) {
                ERROR_EXIT(42, "DEBUG MACRO has not been well configured!");
        }

        backends_functions->init_backend(&devices,
                                         get_nb_dpus_per_run(),
                                         get_dpu_binary(),
                                         &index_seed,
                                         nb_dpu,
                                         ref_genome,
                                         reads_info,
                                         times_ctx,
                                         backends_functions);

        for (int round = 0; round < 3; round++) {
                char filename[1024];
                FILE *fipe1, *fipe2, *fope1, *fope2;
                unsigned int current_pass = 0;
                int8_t *reads_buffer[MAX_NB_PASS];
                int nb_read[MAX_NB_PASS];

                if (DEBUG_ROUND != -1 && DEBUG_ROUND != round) {
                        continue;
                }

                reads_info->delta_neighbour_in_bytes = (SIZE_SEED * round)/4;
                reads_info->size_neighbour_in_32bits_words =
                        (reads_info->size_neighbour_in_bytes-reads_info->delta_neighbour_in_bytes) * 4;

                sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round+1);
                fope1 = fopen(filename, "w");
                sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round+1);
                fope2 = fopen(filename, "w");

                if (round == 0) {
                        sprintf(filename, "%s_PE1.fastq", input_prefix);
                        fipe1 = fopen(filename, "r");
                        sprintf(filename, "%s_PE2.fastq", input_prefix);
                        fipe2 = fopen(filename, "r");
                } else {
                        sprintf(filename, "%s_%d_PE1.fasta", input_prefix, round);
                        fipe1 = fopen(filename, "r");
                        sprintf(filename, "%s_%d_PE2.fasta", input_prefix, round);
                        fipe2 = fopen(filename, "r");
                }

                sprintf(filename, "%s_%d_time.csv", input_prefix, round);
                times_ctx->time_file = fopen(filename, "w");
                fprintf(times_ctx->time_file,
                        "time, write_mram, dispatch_reads, write_reads, compute, "
                        "read_result, map_read, accumulate_read, process_read\n");

                do {
                        assert(current_pass < MAX_NB_PASS);
                        reads_buffer[current_pass] = (int8_t *) malloc(sizeof(int8_t) * MAX_READS_BUFFER * reads_info->size_read);
                        nb_read[current_pass] = get_reads(fipe1, fipe2, reads_buffer[current_pass], times_ctx, reads_info);
                        printf(" - get %d reads\n", nb_read[current_pass] / 2);
                        printf(" - time to get reads      : %7.2lf sec. / %7.2lf sec.\n",
                               times_ctx->get_reads,
                               times_ctx->tot_get_reads);
                } while (nb_read[current_pass++] != 0);

                fclose(fipe1);
                fclose(fipe2);

                map_var_call(round, current_pass - 1,
                             devices, ref_genome, index_seed, dispatch_requests,
                             &variant_list, substitution_list, mapping_coverage,
                             reads_buffer, nb_read,
                             fope1, fope2,
                             reads_info, times_ctx, backends_functions);

                fclose(times_ctx->time_file);
                fclose(fope1);
                fclose(fope2);
        }

        backends_functions->free_backend(devices, nb_dpu);

        create_vcf(input_prefix, ref_genome, &variant_list, substitution_list, mapping_coverage, times_ctx);

        free_variant_tree(variant_list);
        free_genome(ref_genome);
        free_index(index_seed);
        dispatch_free(dispatch_requests, nb_dpu);
        free(substitution_list);
        free(mapping_coverage);
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
        pthread_mutex_init(&times_ctx.time_file_mutex, NULL);

        validate_args(argc, argv);

        printf("%s\n", VERSION);
        print_time();

        reads_info.size_read = get_read_size(get_input_pe1());
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        printf("Information\n");
        printf(" - read size: %d\n", reads_info.size_read);

        setup_dpus_for_target_type(get_target_type());

        if (get_simulation_mode()) {
                backends_functions.init_backend = init_backend_simulation;
                backends_functions.free_backend = free_backend_simulation;
                backends_functions.run_dpu = run_dpu_simulation;
                backends_functions.add_seed_to_requests = add_seed_to_simulation_requests;
                backends_functions.init_vmis = init_vmis_simulation;
                backends_functions.free_vmis = free_vmis_simulation;
                backends_functions.write_vmi = write_vmi_simulation;
                backends_functions.load_mram = load_mram_simulation;
        } else {
                backends_functions.init_backend = init_backend_dpu;
                backends_functions.free_backend = free_backend_dpu;
                backends_functions.run_dpu = run_on_dpu;
                backends_functions.add_seed_to_requests = add_seed_to_dpu_requests;
                backends_functions.init_vmis = init_vmis_dpu;
                backends_functions.free_vmis = free_vmis_dpu;
                backends_functions.write_vmi = write_vmi_dpu;
                backends_functions.load_mram = load_mram_dpu;
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

        pthread_mutex_destroy(&times_ctx.time_file_mutex);
        free_args();

        return 0;
}
