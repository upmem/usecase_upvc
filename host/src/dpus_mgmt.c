/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dispatch.h"
#include "dpu_log.h"
#include "dpus_mgmt.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "upvc.h"

#include "upvc_dpu.h"

#define XSTR(s) STR(s)
#define STR(s) #s

static char *profile;

void setup_dpus_for_target_type(target_type_t target_type)
{
    switch (target_type) {
    case target_type_hw:
        profile = "backend=hw,cycleAccurate=true";
        break;
    default:
        profile = "backend=simulator";
        break;
    }
}

devices_t *dpu_try_alloc_for(unsigned int nb_dpus_per_run, const char *opt_program)
{
    devices_t *devices = (devices_t *)malloc(sizeof(devices_t));
    assert(devices != NULL);

    DPU_ASSERT(dpu_alloc_dpus(profile, nb_dpus_per_run, &devices->ranks, &devices->nb_ranks_per_run, &nb_dpus_per_run));
    set_nb_dpus_per_run(nb_dpus_per_run);

    devices->nb_dpus_per_rank = (unsigned int *)malloc(devices->nb_ranks_per_run * sizeof(unsigned int));
    devices->rank_mram_offset = (unsigned int *)malloc(devices->nb_ranks_per_run * sizeof(unsigned int));
    assert(devices->nb_dpus_per_rank != NULL);
    assert(devices->rank_mram_offset != NULL);

    unsigned int nb_dpus = 0;
    for (unsigned int each_rank = 0; each_rank < devices->nb_ranks_per_run; each_rank++) {
        struct dpu_rank_t *rank = devices->ranks[each_rank];
        DPU_ASSERT(dpu_get_nr_of_dpus_in(rank, &devices->nb_dpus_per_rank[each_rank]));
        DPU_ASSERT(dpu_load_all(rank, opt_program));
        devices->rank_mram_offset[each_rank] = nb_dpus;
        nb_dpus += devices->nb_dpus_per_rank[each_rank];
    }
    printf("%u DPUs allocated\n", nb_dpus);

    struct dpu_t *one_dpu;
    DPU_FOREACH (devices->ranks[0], one_dpu) {
        struct dpu_runtime_context_t *rt_ctx = dpu_get_runtime_context(one_dpu);
        DPU_ASSERT(dpu_get_symbol(rt_ctx, DPU_MRAM_HEAP_POINTER_NAME, &devices->mram_available));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_MRAM_INFO_VAR), &devices->mram_info));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_NB_REQUEST_VAR), &devices->mram_request_info));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_REQUEST_VAR), &devices->mram_requests));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_COMPUTE_TIME_VAR), &devices->mram_compute_time));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_TASKLET_STATS_VAR), &devices->mram_tasklet_stats));
        DPU_ASSERT(dpu_get_symbol(rt_ctx, XSTR(DPU_RESULT_VAR), &devices->mram_result));
        break;
    }

    pthread_mutex_init(&devices->log_mutex, NULL);
    devices->log_file = fopen("upvc_log.txt", "w");

    return devices;
}

void dpu_try_write_mram(unsigned int rank_id, devices_t *devices, uint8_t **mram, size_t *mram_size, int delta_neighbour)
{
    struct dpu_rank_t *rank = devices->ranks[rank_id];
    struct dpu_transfer_mram *matrix, *matrix_delta;
    DPU_ASSERT(dpu_transfer_matrix_allocate(rank, &matrix));
    DPU_ASSERT(dpu_transfer_matrix_allocate(rank, &matrix_delta));

    unsigned int each_dpu = 0;
    struct dpu_t *dpu;
    DPU_FOREACH (rank, dpu) {
        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix, mram[each_dpu], devices->mram_available.size, devices->mram_available.address, 0));
        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix_delta, &delta_neighbour, sizeof(delta_info_t), devices->mram_info.address, 0));
        assert(mram_size[each_dpu] < devices->mram_available.size && "mram is to big to fit in DPU");
        each_dpu++;
    }

    DPU_ASSERT(dpu_copy_to_dpus(rank, matrix));
    DPU_ASSERT(dpu_copy_to_dpus(rank, matrix_delta));
    dpu_transfer_matrix_free(rank, matrix);
    dpu_transfer_matrix_free(rank, matrix_delta);
}

void dpu_try_free(devices_t *devices)
{
    for (unsigned int each_rank = 0; each_rank < devices->nb_ranks_per_run; each_rank++) {
        DPU_ASSERT(dpu_free(devices->ranks[each_rank]));
    }
    pthread_mutex_destroy(&devices->log_mutex);
    fclose(devices->log_file);
    free(devices->ranks);
    free(devices->nb_dpus_per_rank);
    free(devices);
}

void dpu_try_run(unsigned int rank_id, devices_t *devices)
{
    struct dpu_rank_t *rank = devices->ranks[rank_id];
    dpu_api_status_t status = dpu_boot_all(rank, SYNCHRONOUS);
    if (status == DPU_API_DPU_FAULT_ERROR) {
        unsigned int nb_dpus_per_rank;
        uint32_t nb_dpus_running;

        DPU_ASSERT(dpu_get_nr_of_dpus_in(rank, &nb_dpus_per_rank));
        dpu_run_status_t run_status[nb_dpus_per_rank];

        DPU_ASSERT(dpu_get_all_status(rank, run_status, &nb_dpus_running));

        struct dpu_t *dpu;
        unsigned int each_dpu = 0;
        DPU_FOREACH (rank, dpu) {
            switch (run_status[each_dpu]) {
            case DPU_STATUS_IDLE:
            case DPU_STATUS_RUNNING:
                break;
            case DPU_STATUS_ERROR:
                DPU_ASSERT(dpulog_read_for_dpu(dpu, stdout));
                fprintf(stderr, "*** DPU %u.%u reported an error", rank_id, each_dpu);
                break;
            default:
                DPU_ASSERT(dpulog_read_for_dpu(dpu, stdout));
                ERROR_EXIT(11, "*** could not get DPU %u.%u status %u - aborting", rank_id, each_dpu, run_status[each_dpu]);
            }
            each_dpu++;
        }
        exit(10);
    }
    assert(status == DPU_API_SUCCESS && "dpu_boot_all failed");
}

void dpu_try_write_dispatch_into_mram(
    unsigned int rank_id, unsigned int dpu_offset, devices_t *devices, dispatch_request_t *dispatch)
{
    struct dpu_rank_t *rank = devices->ranks[rank_id];
    struct dpu_transfer_mram *matrix_header, *matrix_reads;

    unsigned int nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_of_dpus_in(rank, &nb_dpus_per_rank));
    dispatch_request_t io_header[nb_dpus_per_rank];

    DPU_ASSERT(dpu_transfer_matrix_allocate(rank, &matrix_header));
    DPU_ASSERT(dpu_transfer_matrix_allocate(rank, &matrix_reads));

    dispatch_request_t dummy_dispatch;
    dummy_dispatch.nb_reads = 0;
    dummy_dispatch.reads_area = malloc(sizeof(dpu_request_t) * MAX_DPU_REQUEST);
    assert(dummy_dispatch.reads_area);

    struct dpu_t *dpu;
    unsigned int each_dpu = 0;
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu < nb_dpu) {
            io_header[each_dpu] = dispatch[this_dpu];
        } else {
            io_header[each_dpu] = dummy_dispatch;
        }

        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix_header, &io_header[each_dpu].nb_reads, sizeof(nb_request_t), devices->mram_request_info.address, 0));
        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix_reads, io_header[each_dpu].reads_area, devices->mram_requests.size, devices->mram_requests.address, 0));
        each_dpu++;
    }

    DPU_ASSERT(dpu_copy_to_dpus(rank, matrix_header));
    DPU_ASSERT(dpu_copy_to_dpus(rank, matrix_reads));

    dpu_transfer_matrix_free(rank, matrix_header);
    dpu_transfer_matrix_free(rank, matrix_reads);

    free(dummy_dispatch.reads_area);
}

#define CLOCK_PER_SECONDS (600000000.0)
static void dpu_try_log(unsigned int rank_id, unsigned int dpu_offset, devices_t *devices, struct dpu_transfer_mram *matrix)
{
    struct dpu_rank_t *rank = devices->ranks[rank_id];
    unsigned int nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_of_dpus_in(rank, &nb_dpus_per_rank));
    dpu_compute_time_t compute_time[nb_dpus_per_rank];
    struct dpu_t *dpu;
    unsigned int each_dpu = 0;
    DPU_FOREACH (rank, dpu) {
        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix, &compute_time[each_dpu], sizeof(dpu_compute_time_t), devices->mram_compute_time.address, 0));
        each_dpu++;
    }

    DPU_ASSERT(dpu_copy_from_dpus(rank, matrix));

#ifdef STATS_ON
    /* Collect stats */
    dpu_tasklet_stats_t tasklet_stats[nb_dpus_per_rank][NB_TASKLET_PER_DPU];
    for (unsigned int each_tasklet = 0; each_tasklet < NB_TASKLET_PER_DPU; each_tasklet++) {
        each_dpu = 0;
        DPU_FOREACH (rank, dpu) {
            DPU_ASSERT(
                dpu_transfer_matrix_add_dpu(dpu, matrix, &tasklet_stats[each_dpu][each_tasklet], sizeof(dpu_tasklet_stats_t),
                    (mram_addr_t)(devices->mram_tasklet_stats.address + each_tasklet * sizeof(dpu_tasklet_stats_t)), 0));
            each_dpu++;
        }
        DPU_ASSERT(dpu_copy_from_dpus(rank, matrix));
    }
#endif

    pthread_mutex_lock(&devices->log_mutex);
    fprintf(devices->log_file, "rank %u offset %u\n", rank_id, dpu_offset);
    each_dpu = 0;
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu >= nb_dpu)
            break;
        fprintf(devices->log_file, "LOG DPU=%u TIME=%llu SEC=%.3f\n", this_dpu, (unsigned long long)compute_time[each_dpu],
            (float)compute_time[each_dpu] / CLOCK_PER_SECONDS);
#ifdef STATS_ON
        dpu_tasklet_stats_t agreagated_stats = {
            .nb_reqs = 0,
            .nb_nodp_calls = 0,
            .nb_odpd_calls = 0,
            .nb_results = 0,
            .mram_data_load = 0,
            .mram_result_store = 0,
            .mram_load = 0,
            .mram_store = 0,
            .nodp_time = 0ULL,
            .odpd_time = 0ULL,
        };
        for (unsigned int each_tasklet = 0; each_tasklet < NB_TASKLET_PER_DPU; each_tasklet++) {
            agreagated_stats.nb_reqs += tasklet_stats[each_dpu][each_tasklet].nb_reqs;
            agreagated_stats.nb_nodp_calls += tasklet_stats[each_dpu][each_tasklet].nb_nodp_calls;
            agreagated_stats.nb_odpd_calls += tasklet_stats[each_dpu][each_tasklet].nb_odpd_calls;
            agreagated_stats.nodp_time += (unsigned long long)tasklet_stats[each_dpu][each_tasklet].nodp_time;
            agreagated_stats.odpd_time += (unsigned long long)tasklet_stats[each_dpu][each_tasklet].odpd_time;
            agreagated_stats.nb_results += tasklet_stats[each_dpu][each_tasklet].nb_results;
            agreagated_stats.mram_data_load += tasklet_stats[each_dpu][each_tasklet].mram_data_load;
            agreagated_stats.mram_result_store += tasklet_stats[each_dpu][each_tasklet].mram_result_store;
            agreagated_stats.mram_load += tasklet_stats[each_dpu][each_tasklet].mram_load;
            agreagated_stats.mram_store += tasklet_stats[each_dpu][each_tasklet].mram_store;
        }
        fprintf(devices->log_file, "LOG DPU=%u REQ=%u\n", this_dpu, agreagated_stats.nb_reqs);
        fprintf(devices->log_file, "LOG DPU=%u NODP=%u\n", this_dpu, agreagated_stats.nb_nodp_calls);
        fprintf(devices->log_file, "LOG DPU=%u ODPD=%u\n", this_dpu, agreagated_stats.nb_odpd_calls);
        fprintf(devices->log_file, "LOG DPU=%u NODP_TIME=%llu\n", this_dpu, (unsigned long long)agreagated_stats.nodp_time);
        fprintf(devices->log_file, "LOG DPU=%u ODPD_TIME=%llu\n", this_dpu, (unsigned long long)agreagated_stats.odpd_time);
        fprintf(devices->log_file, "LOG DPU=%u RESULTS=%u\n", this_dpu, agreagated_stats.nb_results);
        fprintf(devices->log_file, "LOG DPU=%u DATA_IN=%u\n", this_dpu, agreagated_stats.mram_data_load);
        fprintf(devices->log_file, "LOG DPU=%u RESULT_OUT=%u\n", this_dpu, agreagated_stats.mram_result_store);
        fprintf(devices->log_file, "LOG DPU=%u LOAD=%u\n", this_dpu, agreagated_stats.mram_load);
        fprintf(devices->log_file, "LOG DPU=%u STORE=%u\n", this_dpu, agreagated_stats.mram_store);

#endif
        DPU_ASSERT(dpulog_read_for_dpu(dpu, devices->log_file));
        each_dpu++;
    }
    fflush(devices->log_file);
    pthread_mutex_unlock(&devices->log_mutex);
}

void dpu_try_get_results_and_log(
    unsigned int rank_id, __attribute__((unused)) unsigned int dpu_offset, devices_t *devices, dpu_result_out_t **result_buffer)
{
    struct dpu_rank_t *rank = devices->ranks[rank_id];
    struct dpu_transfer_mram *matrix;
    DPU_ASSERT(dpu_transfer_matrix_allocate(rank, &matrix));
    uint32_t nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_of_dpus_in(rank, &nb_dpus_per_rank));
    uint8_t *results[nb_dpus_per_rank];
    uint8_t *dummy_results = malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);

    struct dpu_t *dpu;
    unsigned int each_dpu = 0;
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            results[each_dpu] = (uint8_t *)result_buffer[each_dpu];
        } else {
            results[each_dpu] = dummy_results;
        }
        DPU_ASSERT(dpu_transfer_matrix_add_dpu(
            dpu, matrix, results[each_dpu], devices->mram_result.size, devices->mram_result.address, 0));
        each_dpu++;
    }

    DPU_ASSERT(dpu_copy_from_dpus(rank, matrix));

    dpu_try_log(rank_id, dpu_offset, devices, matrix);
    dpu_transfer_matrix_free(rank, matrix);

    free(dummy_results);
}
