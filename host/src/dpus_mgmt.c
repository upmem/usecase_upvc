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

static char *profile = NULL;
static dispatch_request_t dummy_dispatch;
static uint8_t *dummy_results;

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
    struct dpu_program_t *dpu_program;
    devices_t *devices = (devices_t *)malloc(sizeof(devices_t));
    assert(devices != NULL);

    DPU_ASSERT(dpu_alloc(nb_dpus_per_run, profile, &devices->all_ranks));
    DPU_ASSERT(dpu_load(devices->all_ranks, opt_program, &dpu_program));
    DPU_ASSERT(dpu_get_nr_ranks(devices->all_ranks, &devices->nb_ranks_per_run));
    set_nb_dpus_per_run(nb_dpus_per_run);

    devices->nb_dpus_per_rank = (unsigned int *)malloc(devices->nb_ranks_per_run * sizeof(unsigned int));
    devices->rank_mram_offset = (unsigned int *)malloc(devices->nb_ranks_per_run * sizeof(unsigned int));
    devices->ranks = (struct dpu_set_t *)malloc(devices->nb_ranks_per_run * sizeof(struct dpu_set_t));
    assert(devices->nb_dpus_per_rank != NULL);
    assert(devices->rank_mram_offset != NULL);
    assert(devices->ranks != NULL);

    unsigned int nb_dpus = 0;
    unsigned int each_rank;
    struct dpu_set_t rank;
    DPU_RANK_FOREACH (devices->all_ranks, rank, each_rank) {
        devices->ranks[each_rank] = rank;
        DPU_ASSERT(dpu_get_nr_dpus(rank, &devices->nb_dpus_per_rank[each_rank]));
        devices->rank_mram_offset[each_rank] = nb_dpus;
        nb_dpus += devices->nb_dpus_per_rank[each_rank];
    }
    printf("%u DPUs allocated\n", nb_dpus);

    DPU_ASSERT(dpu_get_symbol(dpu_program, DPU_MRAM_HEAP_POINTER_NAME, &devices->mram_available));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_MRAM_INFO_VAR), &devices->mram_info));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_NB_REQUEST_VAR), &devices->mram_request_info));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_REQUEST_VAR), &devices->mram_requests));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_COMPUTE_TIME_VAR), &devices->mram_compute_time));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_TASKLET_STATS_VAR), &devices->mram_tasklet_stats));
    DPU_ASSERT(dpu_get_symbol(dpu_program, XSTR(DPU_RESULT_VAR), &devices->mram_result));

    pthread_mutex_init(&devices->log_mutex, NULL);
    devices->log_file = fopen("upvc_log.txt", "w");

    dummy_dispatch.nb_reads = 0;
    dummy_dispatch.reads_area = malloc(sizeof(dpu_request_t) * MAX_DPU_REQUEST);
    assert(dummy_dispatch.reads_area);
    dummy_results = malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
    assert(dummy_results);

    return devices;
}

void dpu_try_write_mram(unsigned int rank_id, devices_t *devices, uint8_t **mram, size_t *mram_size, int delta_neighbour)
{
    struct dpu_set_t rank = devices->ranks[rank_id];

    unsigned int each_dpu;
    struct dpu_set_t dpu;
    DPU_FOREACH (rank, dpu, each_dpu) {
        assert(mram_size[each_dpu] < devices->mram_available.size && "mram is to big to fit in DPU");
        DPU_ASSERT(dpu_prepare_xfer(dpu, mram[each_dpu]));
    }
    DPU_ASSERT(
        dpu_push_xfer_symbol(rank, DPU_XFER_TO_DPU, devices->mram_available, 0, devices->mram_available.size, DPU_XFER_DEFAULT));

    DPU_ASSERT(dpu_copy_to_symbol(rank, devices->mram_info, 0, &delta_neighbour, devices->mram_info.size));
}

void dpu_try_free(devices_t *devices)
{
    DPU_ASSERT(dpu_free(devices->all_ranks));
    pthread_mutex_destroy(&devices->log_mutex);
    fclose(devices->log_file);
    free(devices->ranks);
    free(devices->nb_dpus_per_rank);
    free(devices);
    free(dummy_dispatch.reads_area);
    free(dummy_results);
}

void dpu_try_run(unsigned int rank_id, devices_t *devices)
{
    struct dpu_set_t rank = devices->ranks[rank_id];
    dpu_error_t status = dpu_launch(rank, DPU_SYNCHRONOUS);
    if (status == DPU_ERR_DPU_FAULT) {
        struct dpu_set_t dpu;
        unsigned int each_dpu;
        DPU_FOREACH (rank, dpu, each_dpu) {
            bool done, fault;
            DPU_ASSERT(dpu_status(dpu, &done, &fault));
            assert(done);
            if (fault) {
                DPU_ASSERT(dpulog_read_for_dpu(dpu.dpu, stdout));
                fprintf(stderr, "*** DPU %u.%u reported an error", rank_id, each_dpu);
            }
        }
    }
    DPU_ASSERT(status);
}

void dpu_try_write_dispatch_into_mram(
    unsigned int rank_id, unsigned int dpu_offset, devices_t *devices, dispatch_request_t *dispatch)
{
    struct dpu_set_t rank = devices->ranks[rank_id];

    unsigned int nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_dpus(rank, &nb_dpus_per_rank));
    dispatch_request_t io_header[nb_dpus_per_rank];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu < nb_dpu) {
            io_header[each_dpu] = dispatch[this_dpu];
        } else {
            io_header[each_dpu] = dummy_dispatch;
        }
    }
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &io_header[each_dpu].nb_reads));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(
        rank, DPU_XFER_TO_DPU, devices->mram_request_info, 0, devices->mram_request_info.size, DPU_XFER_DEFAULT));

    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, io_header[each_dpu].reads_area));
    }
    DPU_ASSERT(
        dpu_push_xfer_symbol(rank, DPU_XFER_TO_DPU, devices->mram_requests, 0, devices->mram_requests.size, DPU_XFER_DEFAULT));
}

#define CLOCK_PER_SECONDS (600000000.0)
static void dpu_try_log(unsigned int rank_id, unsigned int dpu_offset, devices_t *devices)
{
    struct dpu_set_t rank = devices->ranks[rank_id];
    unsigned int nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_dpus(rank, &nb_dpus_per_rank));
    dpu_compute_time_t compute_time[nb_dpus_per_rank];
    struct dpu_set_t dpu;
    unsigned int each_dpu;
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &compute_time[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(
        rank, DPU_XFER_FROM_DPU, devices->mram_compute_time, 0, sizeof(dpu_compute_time_t), DPU_XFER_DEFAULT));

#ifdef STATS_ON
    /* Collect stats */
    dpu_tasklet_stats_t tasklet_stats[nb_dpus_per_rank][NR_TASKLETS];
    for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
        DPU_FOREACH (rank, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &tasklet_stats[each_dpu][each_tasklet]));
        }
        DPU_ASSERT(dpu_push_xfer_symbol(dpu, DPU_XFER_FROM_DPU, devices->mram_tasklet_stats,
            each_tasklet * sizeof(dpu_tasklet_stats_t), sizeof(dpu_tasklet_stats_t), DPU_XFER_DEFAULT));
    }
#endif

    pthread_mutex_lock(&devices->log_mutex);
    fprintf(devices->log_file, "rank %u offset %u\n", rank_id, dpu_offset);
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu, each_dpu) {
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
        for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
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
        DPU_ASSERT(dpu_log_read(dpu, devices->log_file));
    }
    fflush(devices->log_file);
    pthread_mutex_unlock(&devices->log_mutex);
}

void dpu_try_get_results_and_log(
    unsigned int rank_id, __attribute__((unused)) unsigned int dpu_offset, devices_t *devices, dpu_result_out_t **result_buffer)
{
    struct dpu_set_t rank = devices->ranks[rank_id];
    uint32_t nb_dpus_per_rank;
    DPU_ASSERT(dpu_get_nr_dpus(rank, &nb_dpus_per_rank));
    uint8_t *results[nb_dpus_per_rank];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            results[each_dpu] = (uint8_t *)result_buffer[each_dpu];
        } else {
            results[each_dpu] = dummy_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results[each_dpu]));
    }
    DPU_ASSERT(
        dpu_push_xfer_symbol(rank, DPU_XFER_FROM_DPU, devices->mram_result, 0, devices->mram_result.size, DPU_XFER_DEFAULT));

    dpu_try_log(rank_id, dpu_offset, devices);
}
