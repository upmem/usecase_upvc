/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "parse_args.h"
#include "dpus.h"
#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "upvc_dpu.h"
#include "processread.h"
#include "vartree.h"
#include "vcf.h"
#include "genome.h"
#include "upvc.h"
#include "mdpu.h"

static void run_dpu_simulation(unsigned int nb_dpu,
                               times_ctx_t *times_ctx,
                               reads_info_t *reads_info)
{
        unsigned int numdpu;
        double t1, t2;
        pthread_t thread_id[nb_dpu];
        align_on_dpu_arg_t thread_args[nb_dpu];

        t1 = my_clock();

        for (numdpu = 0; numdpu < nb_dpu; numdpu++) {
                thread_args[numdpu].reads_info = *reads_info;
                thread_args[numdpu].numdpu = numdpu;
        }

        for (numdpu = 0; numdpu < nb_dpu; numdpu++) {
                pthread_create(&thread_id[numdpu], NULL, align_on_dpu, &thread_args[numdpu]);
        }
        for (numdpu = 0; numdpu < nb_dpu; numdpu++) {
                pthread_join(thread_id[numdpu], NULL);
        }

        t2 = my_clock();
        times_ctx->map_read = t2 - t1;
        times_ctx->tot_map_read += t2 - t1;
}

static void run_on_dpu(dispatch_t dispatch,
                       const char *dpu_binary,
                       unsigned int nb_dpu,
                       times_ctx_t *times_ctx,
                       reads_info_t *reads_info)
{
        double t1, t2;
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
        uint64_t t0s[nb_dpus_per_run];

        t1 = my_clock();

        for (unsigned int first_dpu = 0; first_dpu < nb_dpu; first_dpu += nb_dpus_per_run) {
                devices_t devices = dpu_try_alloc_for(nb_dpus_per_run, dpu_binary);
                if (devices == NULL) {
                        ERROR("Unable to alloc devices!");
                }

                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                printf("() write MRAM #%d %u reads\n", this_dpu, dispatch[this_dpu].nb_reads);
                                dpu_try_load_mram_number(this_dpu, each_dpu, devices, reads_info);
                                dpu_try_write_dispatch_into_mram(each_dpu, devices,
                                                                 dispatch[this_dpu].nb_reads,
                                                                 dispatch[this_dpu].reads_area,
                                                                 reads_info);
                                /* DEBUG */
                                /* if (each_dpu == 0) dpu_try_backup_mram(each_dpu, devices, "mram_dbg.bin"); */
                        }
                }

                unsigned int nb_booted_dpus = 0;
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
#ifdef TEST_BACKUP_MRAM_ON_DPU_NR
                                if (this_dpu == TEST_BACKUP_MRAM_ON_DPU_NR) {
                                        dpu_try_backup_mram(this_dpu, dpus[each_dpu], "mram.bck");
                                }
#endif
                                printf("() boot DPU #%d\n", this_dpu);
                                t0s[each_dpu] = dpu_try_run(each_dpu, devices);
                                nb_booted_dpus++;
                        }
                }

                /* Wait for DPUs to complete: use a bitfield to mark every DPU stopped. */
                uint32_t cmask = 0, mask_ok = 0;
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                mask_ok |= (1 << each_dpu);
                        }
                }
                /* int debug_count = 0, debug_iter_count = 0; */
                do {
                        /* debug_count++; */
                        /* if (debug_count == (512 << 20)) { */
                        /*         debug_iter_count++; */
                        /*         if (debug_iter_count == 10) { */
                        /*                 fprintf(stderr, "no answer from DPU! aborting simulation\n"); */
                        /*                 exit(22); */
                        /*         } */
                        /*         debug_count = 0; */
                        /* } */

                        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                                uint32_t mask = (uint32_t) (1 << each_dpu);
                                unsigned int this_dpu = first_dpu + each_dpu;
                                if (!(cmask & mask)) {
                                        if (dispatch[this_dpu].nb_reads != 0) {
                                                if (dpu_try_check_status(each_dpu, devices)) {
                                                        cmask |= mask;
                                                        printf("DPU #%u completed\n", this_dpu);
                                                }
                                        }
                                }
                        }
                } while (mask_ok != cmask);

                /* Gather results and free DPUs */
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                dpu_try_log(each_dpu, devices, t0s[each_dpu]);
                                dpu_result_out_t *results = dpu_try_get_results(this_dpu, devices);
                                int i;
                                for (i = 0; results[i].num != -1; i++) {
                                        if (i == MAX_DPU_RESULTS) {
                                                ERROR_EXIT(22, "BUG! no EOR marker found when parsing DPU results!\n");
                                        }

                                        long coords = ((long) results[i].seed_nr) |
                                                (((long) (results[i].seq_nr)) << 32);
                                        write_result(this_dpu, i, results[i].num, coords, results[i].score);
#ifdef TEST_PRINT_FIRST_DPU_RESULTS
                                        if (each_dpu == 0) {
                                                printf("num=%u seed=%u seq=%u score=%u\n",
                                                       results[i].num,
                                                       results[i].seed_nr,
                                                       results[i].seq_nr,
                                                       results[i].score);
                                        }
#endif /* TEST_PRINT_FIRST_DPU_RESULTS */
                                }
                                write_result(this_dpu, i, (unsigned int) -1, -1L, (unsigned int) -1);

                                free(results);
                        }
                }

                dpu_try_free(devices);
        }

        t2 = my_clock();
        times_ctx->map_read = t2 - t1;
        times_ctx->tot_map_read += t2 - t1;
}

static void run_dpu(dispatch_t dispatch,
                    char *dpu_binary,
                    unsigned int nb_dpu,
                    times_ctx_t *times_ctx,
                    reads_info_t *reads_info,
                    bool simulation_mode)
{
        if (simulation_mode) {
                run_dpu_simulation(nb_dpu, times_ctx, reads_info);
        } else {
                run_on_dpu(dispatch, dpu_binary, nb_dpu, times_ctx, reads_info);
        }
}

static int map_var_call(char *filename_prefix,
                        int round,
                        genome_t *ref_genome,
                        index_seed_t **index_seed,
                        variant_tree_t **variant_list,
                        int *substitution_list,
                        int8_t *mapping_coverage,
                        reads_info_t *reads_info,
                        times_ctx_t *times_ctx,
                        bool simulation_mode)
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

                dispatch = dispatch_read(index_seed,
                                         reads_buffer,
                                         nb_read,
                                         nb_dpu,
                                         times_ctx,
                                         reads_info,
                                         simulation_mode);
                printf(" - time to dispatch reads : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->dispatch_read,
                       times_ctx->tot_dispatch_read);

                run_dpu(dispatch, get_dpu_binary(), nb_dpu, times_ctx, reads_info, simulation_mode);
                printf(" - time to map reads      : %7.2lf sec. / %7.2lf sec.\n",
                       times_ctx->map_read,
                       times_ctx->tot_map_read);

                dispatch_free(dispatch, nb_dpu);

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

static index_seed_t **reload_mram_images_and_seeds(reads_info_t *reads_info)
{
        index_seed_t **index_seed;
        mram_info_t *mram = mram_create(reads_info);
        unsigned int nb_dpus = get_nb_dpu();
        /* Will unwrap the MRAM contents into the MDPU's PMEM. */

        for (unsigned int each_dpu = 0; each_dpu < nb_dpus; each_dpu++) {
                mram_reset(mram, reads_info);
                mram_load(mram, each_dpu);
                malloc_neighbour_idx(each_dpu, mram->nb_nbr, reads_info);
                write_neighbours_and_coordinates(each_dpu, mram->nb_nbr, mram_neighbours_area(mram), reads_info);
        }

        index_seed = load_index_seeds();

        mram_free(mram);

        return index_seed;
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

static void load_index_save_genome(reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);
        index_seed_t **index_seed = index_genome(ref_genome,
                                                 get_nb_dpu(),
                                                 times_ctx,
                                                 reads_info,
                                                 get_simulation_mode());
        save_index_seeds(index_seed);

        free_genome(ref_genome);
        free_index(index_seed);
}

static void do_mapping(bool simulation_mode, reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        char *input_prefix = get_input_path();
        unsigned int nb_dpu = get_nb_dpu();
        variant_tree_t *variant_list = NULL;
        genome_t *ref_genome = get_genome(get_input_fasta(), times_ctx);

        int8_t *mapping_coverage = (int8_t *) calloc(sizeof(int8_t), ref_genome->fasta_file_size);
        int *substitution_list = (int *) calloc(sizeof(int), ref_genome->fasta_file_size);

        index_seed_t **index_seed;
        if (simulation_mode) {
                index_seed = index_genome(ref_genome, nb_dpu, times_ctx, reads_info, simulation_mode);
        } else {
                malloc_dpu(reads_info, nb_dpu);
                index_seed = reload_mram_images_and_seeds(reads_info);
        }

        map_var_call(input_prefix, 0, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     reads_info, times_ctx, simulation_mode);
        map_var_call(input_prefix, 1, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     reads_info, times_ctx, simulation_mode);
        map_var_call(input_prefix, 2, ref_genome, index_seed,
                     &variant_list, substitution_list, mapping_coverage,
                     reads_info, times_ctx, simulation_mode);

        create_vcf(input_prefix, ref_genome, &variant_list, substitution_list, mapping_coverage, times_ctx);

        free_variant_tree(variant_list);
        free_genome(ref_genome);
        free_index(index_seed);
        free(substitution_list);
        free(mapping_coverage);
        free_dpu(nb_dpu);
}

int main(int argc, char *argv[])
{
        reads_info_t reads_info;
        times_ctx_t times_ctx;

        memset(&times_ctx, 0, sizeof(times_ctx_t));

        validate_args(argc, argv);

        printf("%s\n", VERSION);

        reads_info.size_read = get_read_size(get_input_pe1());
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        printf("Information\n");
        printf(" - read size: %d\n", reads_info.size_read);

        setup_dpus_for_target_type(get_target_type());

        switch(get_goal()) {
        case goal_index:
                load_index_save_genome(&reads_info, &times_ctx);
                break;
        case goal_check:
                reload_and_verify_mram_images(&reads_info);
                break;
        case goal_map:
                do_mapping(get_simulation_mode(), &reads_info, &times_ctx);
                break;
        case goal_unknown:
        default:
                ERROR_EXIT(23, "goal has not been specified!");
        }

        free_args();

        return 0;
}
