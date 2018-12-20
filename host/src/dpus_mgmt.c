/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "upvc.h"
#include "mram_dpu.h"
#include "dpus_mgmt.h"
#include "dpulog.h"
#include "dispatch.h"
#include "parse_args.h"

#define DPU_SLICE_SIZE (8)

/* #define LOG_DPUS */
#ifdef LOG_DPUS
static dpu_logging_config_t logging_config = {
                                              .source = KTRACE,
                                              .destination_directory_name = "."
};
#define p_logging_config &logging_config
#define log_dpu(dpu, out) dpulog_read_for_dpu(dpu, out)
#else
#define p_logging_config NULL
#define log_dpu(dpu, log) do {} while(0)
#endif /* LOG_DPUS */

static struct dpu_param param = {
                                 .type = FUNCTIONAL_SIMULATOR,
                                 .profile = "",
                                 .on_boot = NULL,
                                 .logging_config = p_logging_config
};

void setup_dpus_for_target_type(target_type_t target_type)
{
        switch (target_type) {
        case target_type_fpga:
                param.type = HW;
                break;
        default:
                param.type = FUNCTIONAL_SIMULATOR;
                break;
        }
}

devices_t *dpu_try_alloc_for(unsigned int nb_dpus_per_run, const char *opt_program)
{
        dpu_api_status_t status;
        devices_t *devices = (devices_t *) malloc(sizeof(devices_t));
        assert(devices != NULL);

        devices->nb_dpus = nb_dpus_per_run;
        status = dpu_get_nr_of_dpus_for(&param, &(devices->nb_dpus_per_rank));
        assert(status == DPU_API_SUCCESS && "dpu_get_nr_of_dpus_for failed");

        devices->dpus = (dpu_t *) calloc(nb_dpus_per_run, sizeof(dpu_t));
        assert(devices->dpus != NULL);

        devices->mram_info = (mram_info_t *) malloc(get_nb_dpu() * sizeof(mram_info_t));
        assert(devices->mram_info != NULL);

        if (nb_dpus_per_run % devices->nb_dpus_per_rank != 0) {
                ERROR_EXIT(5, "*** number of DPUs per run is not a multiple of the DPUs in a rank - aborting");
        }
        devices->nb_ranks_per_run = nb_dpus_per_run / devices->nb_dpus_per_rank;
        devices->ranks = (dpu_rank_t *) calloc(devices->nb_ranks_per_run, sizeof(dpu_rank_t));
        assert(devices->ranks != NULL);

        for (unsigned int each_rank = 0, each_dpu = 0; each_dpu < nb_dpus_per_run; each_rank++) {
                status = dpu_alloc(&param, &(devices->ranks[each_rank]));
                assert(status == DPU_API_SUCCESS && "dpu_alloc failed");
                for (unsigned int each_member = 0;
                     (each_member < devices->nb_dpus_per_rank) && (each_dpu < nb_dpus_per_run);
                     each_member++, each_dpu++) {
                        devices->dpus[each_dpu] = dpu_get_id(devices->ranks[each_rank], each_member);
                }
        }
        for (unsigned int each_dpu = 0; each_dpu < get_nb_dpu(); each_dpu++) {
                mram_load_info(&devices->mram_info[each_dpu], each_dpu);
        }

        for (unsigned int each_rank = 0; each_rank < devices->nb_ranks_per_run; each_rank++) {
                status = dpu_load_all(devices->ranks[each_rank], opt_program);
                assert(status == DPU_API_SUCCESS && "dpu_load_all failed");
        }

        return devices;
}

void dpu_try_write_mram(unsigned int rank_id, devices_t *devices, mram_info_t **mram)
{
        dpu_api_status_t status;
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;
        dpu_rank_t rank = devices->ranks[rank_id];
        dpu_transfer_mram_t *matrix;
        status = dpu_transfer_matrix_allocate(rank, &matrix);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_allocate failed");
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                status = dpu_transfer_matrix_add_dpu(rank,
                                                     matrix,
                                                     each_dpu/DPU_SLICE_SIZE,
                                                     each_dpu%DPU_SLICE_SIZE,
                                                     mram[each_dpu],
                                                     sizeof(mram_info_t) + mram[each_dpu]->total_nbr_size,
                                                     MRAM_INFO_ADDR,
                                                     0);
                assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
        }

        status = dpu_copy_to_dpus(rank, matrix);
        assert(status == DPU_API_SUCCESS && "dpu_copy_to_dpus failed");
        status = dpu_transfer_matrix_free(rank, matrix);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_free failed");
}

void dpu_try_free(devices_t *devices)
{
        for (unsigned int each_rank = 0; each_rank < devices->nb_ranks_per_run; each_rank++) {
                dpu_free(devices->ranks[each_rank]);
        }
        free(devices->ranks);
        free(devices->dpus);
        free(devices->mram_info);
        free(devices);
}

void dpu_try_run(unsigned int rank_id, devices_t *devices)
{
        dpu_api_status_t status = dpu_boot_all(devices->ranks[rank_id], ASYNCHRONOUS);
        assert(status == DPU_API_SUCCESS && "dpu_boot_all failed");
}

bool dpu_try_check_status(unsigned int rank_id, devices_t *devices)
{
        dpu_api_status_t status;
        dpu_run_status_t run_status[devices->nb_dpus_per_rank];
        uint32_t nb_dpus_running;
        status = dpu_get_all_status(devices->ranks[rank_id], run_status, &nb_dpus_running);
        assert(status == DPU_API_SUCCESS && "dpu_get_all_status failed");

        for (unsigned int each_dpu = 0; each_dpu < devices->nb_dpus_per_rank; each_dpu++) {
                switch (run_status[each_dpu]) {
                case DPU_STATUS_IDLE:
                case DPU_STATUS_RUNNING:
                        continue;
                case DPU_STATUS_ERROR:
                        log_dpu(devices->dpus[each_dpu], stdout);
                        ERROR_EXIT(10, "*** DPU %u reported an error - aborting", each_dpu);
                default:
                        log_dpu(devices->dpus[each_dpu], stdout);
                        ERROR_EXIT(11, "*** could not get DPU status %u - aborting", each_dpu);
                }
        }
        return (nb_dpus_running == 0);
}

void dpu_try_write_dispatch_into_mram(unsigned int rank_id,
                                      unsigned int dpu_offset,
                                      devices_t *devices,
                                      dispatch_request_t *dispatch,
                                      reads_info_t *reads_info)
{
        dpu_api_status_t status;
        dpu_rank_t rank = devices->ranks[rank_id];
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;
        dpu_transfer_mram_t *matrix_header, *matrix_reads;
        request_info_t io_header[nb_dpus_per_rank];
        unsigned int io_len[nb_dpus_per_rank];

        status = dpu_transfer_matrix_allocate(rank, &matrix_header);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_allocate failed");
        status = dpu_transfer_matrix_allocate(rank, &matrix_reads);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_allocate failed");

        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                mram_info_t *mram = &devices->mram_info[each_dpu + dpu_offset];
                unsigned int nb_reads = dispatch[each_dpu + dpu_offset].nb_reads;
                io_len[each_dpu] = nb_reads * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes);
                io_header[each_dpu].nb_reads = nb_reads;
                io_header[each_dpu].magic = 0xcdefabcd;

                if ((DPU_REQUEST_ADDR(mram) - DPU_INPUTS_ADDR + io_len[each_dpu]) > DPU_INPUTS_SIZE) {
                        ERROR_EXIT(12,
                                   "*** will exceed MRAM limit if writing reads on DPU number %u - aborting!",
                                   each_dpu + dpu_offset);
                }
                status = dpu_transfer_matrix_add_dpu(rank,
                                                     matrix_header,
                                                     each_dpu/DPU_SLICE_SIZE,
                                                     each_dpu%DPU_SLICE_SIZE,
                                                     &io_header[each_dpu],
                                                     sizeof(request_info_t),
                                                     (mram_addr_t) DPU_REQUEST_INFO_ADDR(mram),
                                                     0);
                assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
                status = dpu_transfer_matrix_add_dpu(rank,
                                                     matrix_reads,
                                                     each_dpu/DPU_SLICE_SIZE,
                                                     each_dpu%DPU_SLICE_SIZE,
                                                     dispatch[each_dpu + dpu_offset].reads_area,
                                                     io_len[each_dpu],
                                                     (mram_addr_t) DPU_REQUEST_ADDR(mram),
                                                     0);
                assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
        }

        status = dpu_copy_to_dpus(rank, matrix_header);
        assert(status == DPU_API_SUCCESS && "dpu_copy_to_dpus failed");
        status = dpu_copy_to_dpus(rank, matrix_reads);
        assert(status == DPU_API_SUCCESS && "dpu_copy_to_dpus failed");

        status = dpu_transfer_matrix_free(rank, matrix_header);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_free failed");
        status = dpu_transfer_matrix_free(rank, matrix_reads);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_free failed");
}

#define CLOCK_PER_SEC (600000000.0)
static void dpu_try_log(unsigned int rank_id,
                        unsigned int dpu_offset,
                        devices_t *devices,
                        dpu_transfer_mram_t *matrix)
{
        dpu_api_status_t status;
        dpu_rank_t rank = devices->ranks[rank_id];
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;
        dpu_compute_time_t compute_time[nb_dpus_per_rank];
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                status = dpu_transfer_matrix_add_dpu(rank,
                                                     matrix,
                                                     each_dpu/DPU_SLICE_SIZE,
                                                     each_dpu%DPU_SLICE_SIZE,
                                                     &compute_time[each_dpu],
                                                     sizeof(dpu_compute_time_t),
                                                     (mram_addr_t) DPU_COMPUTE_TIME_ADDR,
                                                     0);
                assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
        }

        status = dpu_copy_from_dpus(rank, matrix);
        assert(status == DPU_API_SUCCESS && "dpu_copy_from_dpus failed");

#ifdef STATS_ON
        /* Collect stats */
        dpu_tasklet_stats_t tasklet_stats[nb_dpus_per_rank][NB_TASKLET_PER_DPU];
        for (unsigned int each_tasklet = 0; each_tasklet < NB_TASKLET_PER_DPU; each_tasklet++) {
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                        status = dpu_transfer_matrix_add_dpu(rank,
                                                             matrix,
                                                             each_dpu/DPU_SLICE_SIZE,
                                                             each_dpu%DPU_SLICE_SIZE,
                                                             &tasklet_stats[each_dpu][each_tasklet],
                                                             sizeof(dpu_tasklet_stats_t),
                                                             (mram_addr_t)
                                                             (DPU_TASKLET_STATS_ADDR + each_tasklet * sizeof(dpu_tasklet_stats_t)),
                                                             0);
                        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
                }
                status = dpu_copy_from_dpus(rank, matrix);
                assert(status == DPU_API_SUCCESS && "dpu_copy_from_dpus failed");
        }
#endif

        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                unsigned int this_dpu = each_dpu + dpu_offset + rank_id * nb_dpus_per_rank;
                printf("LOG DPU=%u TIME=%llu SEC=%.3f\n",
                       this_dpu,
                       (unsigned long long)compute_time[each_dpu],
                       (float)compute_time[each_dpu] / CLOCK_PER_SEC);
#ifdef STATS_ON
                for (unsigned int each_tasklet = 0; each_tasklet < NB_TASKLET_PER_DPU; each_tasklet++) {
                        /* TODO: show logs */
                        printf("LOG DPU=%u TID=%u REQ=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].nb_reqs);
                        printf("LOG DPU=%u TID=%u NODP=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].nb_nodp_calls);
                        printf("LOG DPU=%u TID=%u ODPD=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].nb_odpd_calls);
                        printf("LOG DPU=%u TID=%u RESULTS=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].nb_results);
                        printf("LOG DPU=%u TID=%u DATA_IN=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].mram_data_load);
                        printf("LOG DPU=%u TID=%u RESULT_OUT=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].mram_result_store);
                        printf("LOG DPU=%u TID=%u LOAD=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].mram_load);
                        printf("LOG DPU=%u TID=%u STORE=%u\n",
                               this_dpu, each_tasklet, tasklet_stats[each_dpu][each_tasklet].mram_store);
                }
#endif
                log_dpu(devices->dpus[this_dpu], stdout);
        }
        fflush(stdout);
}

void dpu_try_get_results_and_log(unsigned int rank_id, unsigned int dpu_offset, devices_t *devices, dpu_result_out_t **result_buffer)
{
        dpu_api_status_t status;
        dpu_rank_t rank = devices->ranks[rank_id];
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;
        dpu_transfer_mram_t *matrix;
        status = dpu_transfer_matrix_allocate(rank, &matrix);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_allocate failed");
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                status = dpu_transfer_matrix_add_dpu(rank,
                                                     matrix,
                                                     each_dpu/DPU_SLICE_SIZE,
                                                     each_dpu%DPU_SLICE_SIZE,
                                                     result_buffer[each_dpu],
                                                     DPU_RESULT_SIZE,
                                                     (mram_addr_t) DPU_RESULT_ADDR,
                                                     0);
                assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_add_dpu failed");
        }
        status = dpu_copy_from_dpus(rank, matrix);
        assert(status == DPU_API_SUCCESS && "dpu_copy_from_dpus failed");

        dpu_try_log(rank_id, dpu_offset, devices, matrix);

        status = dpu_transfer_matrix_free(rank, matrix);
        assert(status == DPU_API_SUCCESS && "dpu_transfer_matrix_free failed");
}

void dpu_try_backup_mram(unsigned int tid, devices_t *devices, const char *file_name)
{
        FILE *f = fopen(file_name, "wb");
        printf("saving DPU %u MRAM into '%s'\n", tid, file_name);
        if (f != NULL) {
                uint8_t *mram = (uint8_t *) malloc(MRAM_SIZE);
                dpu_api_status_t status = dpu_copy_from_individual(devices->dpus[tid], 0, mram, MRAM_SIZE);
                assert(status == DPU_API_SUCCESS && "dpu_copy_from_individual failed");
                fwrite(mram, 1, MRAM_SIZE, f);
                free(mram);
                fclose(f);
        } else {
                WARNING("failed to backup DPU %u MRAM into '%s'\n", tid, file_name);
        }
}
