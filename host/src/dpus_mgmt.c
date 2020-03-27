/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <dpu.h>
#include <dpu_log.h>

#include "dispatch.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "upvc.h"
#include "upvc_dpu.h"

uint32_t compute_checksum(uint32_t *buffer, size_t size_in_bytes) {
    uint32_t checksum = 0;
    for (unsigned int i = 0; i < (size_in_bytes / sizeof(uint32_t)); i++) {
        checksum += buffer[i];
    }
    return checksum;
}

#define MAX(x, y) ((x) > (y) ? (x) : (y))

DPU_INCBIN(upvc_dpu_program, DPU_BINARY);

typedef struct {
    unsigned int nb_ranks_per_run;
    unsigned int nb_dpus_per_rank[NB_RANKS_MAX];
    unsigned int rank_mram_offset[NB_RANKS_MAX];
    struct dpu_set_t all_ranks;
    struct dpu_set_t ranks[NB_RANKS_MAX];
    pthread_mutex_t log_mutex;
    FILE *log_file;
} devices_t;

static devices_t devices = { .nb_ranks_per_run = 1 };

#define CHECKSUM_REF_ELEM (640*250*3)
static uint32_t *checksums_ref;
extern volatile int dpus_mgmt_round;
extern volatile int dpus_mgmt_pass[10];
#if 1
#define GENERATE_CHECKUM_REF
#else
#define ASSERT_CHECKSUM_REF
#endif

unsigned int dpu_get_nb_ranks_per_run() { return devices.nb_ranks_per_run; }

unsigned int dpu_try_alloc_for()
{
    unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
    const char *profile = "cycleAccurate=true";

    DPU_ASSERT(dpu_alloc(nb_dpus_per_run, profile, &devices.all_ranks));
    DPU_ASSERT(dpu_load_from_incbin(devices.all_ranks, &upvc_dpu_program, NULL));
    DPU_ASSERT(dpu_get_nr_ranks(devices.all_ranks, &devices.nb_ranks_per_run));
    DPU_ASSERT(dpu_get_nr_dpus(devices.all_ranks, &nb_dpus_per_run));
    assert(devices.nb_ranks_per_run <= NB_RANKS_MAX);
    set_nb_dpus_per_run(nb_dpus_per_run);

    unsigned int nb_dpus = 0;
    unsigned int each_rank;
    struct dpu_set_t rank;
    DPU_RANK_FOREACH (devices.all_ranks, rank, each_rank) {
        devices.ranks[each_rank] = rank;
        DPU_ASSERT(dpu_get_nr_dpus(rank, &devices.nb_dpus_per_rank[each_rank]));
        devices.rank_mram_offset[each_rank] = nb_dpus;
        nb_dpus += devices.nb_dpus_per_rank[each_rank];
    }
    assert(nb_dpus == nb_dpus_per_run);
    printf("%u DPUs allocated\n", nb_dpus);

    pthread_mutex_init(&devices.log_mutex, NULL);
    devices.log_file = fopen("upvc_log.txt", "w");

    checksums_ref = malloc(sizeof(uint32_t) * CHECKSUM_REF_ELEM);
    assert(checksums_ref != NULL);
#ifdef ASSERT_CHECKSUM_REF
    FILE *checksum_file = fopen("checksum_ref.bin", "r");
    fread(checksums_ref, sizeof(uint32_t), CHECKSUM_REF_ELEM, checksum_file);
    fclose(checksum_file);
#endif

    return devices.nb_ranks_per_run;
}

void dpu_try_write_mram(unsigned int rank_id, unsigned int dpu_offset, delta_info_t delta_neighbour)
{
    unsigned int nb_dpu = get_nb_dpu();
    struct dpu_set_t rank = devices.ranks[rank_id];
    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    uint8_t *mram[nb_dpus_per_rank];
    uint32_t mram_size[nb_dpus_per_rank];

    memset(mram_size, 0, sizeof(uint32_t) * nb_dpus_per_rank);

    unsigned int each_dpu;
    struct dpu_set_t dpu;
    dpu_offset += devices.rank_mram_offset[rank_id];
    DPU_FOREACH (rank, dpu, each_dpu) {
        mram[each_dpu] = (uint8_t *)malloc(MRAM_SIZE);
        assert(mram[each_dpu] != NULL);
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            mram_size[each_dpu] = mram_load(mram[each_dpu], this_dpu);
        }
    }

    uint32_t checksums[nb_dpus_per_rank];
    DPU_FOREACH (rank, dpu, each_dpu) {
        checksums[each_dpu] = compute_checksum((uint32_t *)(mram[each_dpu]), mram_size[each_dpu]);
        DPU_ASSERT(dpu_prepare_xfer(dpu, &checksums[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, "mram_checksum",  0, sizeof(uint32_t), DPU_XFER_DEFAULT));
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &mram_size[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, "mram_size", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));

    unsigned int max_mram_size = 0;
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, mram[each_dpu]));
        max_mram_size = MAX(max_mram_size, mram_size[each_dpu]);
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, 0, max_mram_size, DPU_XFER_DEFAULT));

    DPU_ASSERT(dpu_copy_to(rank, XSTR(DPU_MRAM_INFO_VAR), 0, &delta_neighbour, sizeof(delta_neighbour)));

    DPU_FOREACH (rank, dpu, each_dpu) {
        free(mram[each_dpu]);
    }
}

void dpu_try_free()
{
    DPU_ASSERT(dpu_free(devices.all_ranks));
    pthread_mutex_destroy(&devices.log_mutex);
    fclose(devices.log_file);
#ifdef GENERATE_CHECKSUM_REF
    FILE *checksum_file = fopen("checksum_ref.bin", "w");
    fwrite(checksums_ref, sizeof(uint32_t), CHECKSUM_REF_ELEM, checksum_file);
    fclose(checksum_file);
#endif
    free(checksums_ref);
}

void dpu_try_run(unsigned int rank_id)
{
    struct dpu_set_t rank = devices.ranks[rank_id];
    dpu_error_t status = dpu_launch(rank, DPU_SYNCHRONOUS);
    if (status == DPU_ERR_DPU_FAULT) {
        struct dpu_set_t dpu;
        unsigned int each_dpu;
        DPU_FOREACH (rank, dpu, each_dpu) {
            bool done, fault;
            DPU_ASSERT(dpu_status(dpu, &done, &fault));
            assert(done);
            if (fault) {
                DPU_ASSERT(dpu_log_read(dpu, stdout));
                fprintf(stderr, "*** DPU %u.%u reported an error", rank_id, each_dpu);
            }
        }
    }
    DPU_ASSERT(status);
}

void dpu_try_write_dispatch_into_mram(unsigned int rank_id, unsigned int dpu_offset)
{
    static const dpu_request_t dummy_reads_area[MAX_DPU_REQUEST];
    static dispatch_request_t dummy_dispatch = { .nb_reads = 0, .reads_area = (int8_t *)dummy_reads_area };
    struct dpu_set_t rank = devices.ranks[rank_id];

    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    dispatch_request_t *io_header[nb_dpus_per_rank];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = get_nb_dpu();
    unsigned int max_dispatch_size = 0;
    dpu_offset += devices.rank_mram_offset[rank_id];
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu < nb_dpu) {
            io_header[each_dpu] = dispatch_get(this_dpu);
        } else {
            io_header[each_dpu] = &dummy_dispatch;
        }
        max_dispatch_size = MAX(max_dispatch_size, io_header[each_dpu]->nb_reads * sizeof(dpu_request_t));
    }

    uint32_t checksums[nb_dpus_per_rank];
    DPU_FOREACH(rank, dpu, each_dpu) {
        checksums[each_dpu] = compute_checksum(
            (uint32_t *)(io_header[each_dpu]->reads_area), io_header[each_dpu]->nb_reads * sizeof(dpu_request_t));
        DPU_ASSERT(dpu_prepare_xfer(dpu, &checksums[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, "input_checksum", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));

    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &io_header[each_dpu]->nb_reads));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, XSTR(DPU_NB_REQUEST_VAR), 0, sizeof(nb_request_t), DPU_XFER_DEFAULT));

    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, io_header[each_dpu]->reads_area));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, XSTR(DPU_REQUEST_VAR), 0, max_dispatch_size, DPU_XFER_DEFAULT));
}

static void dpu_try_log(unsigned int rank_id, unsigned int dpu_offset)
{
    struct dpu_set_t rank = devices.ranks[rank_id];
    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    dpu_compute_time_t compute_time[nb_dpus_per_rank];
    struct dpu_set_t dpu;
    unsigned int each_dpu;
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &compute_time[each_dpu]));
    }
    DPU_ASSERT(
        dpu_push_xfer(rank, DPU_XFER_FROM_DPU, XSTR(DPU_COMPUTE_TIME_VAR), 0, sizeof(dpu_compute_time_t), DPU_XFER_DEFAULT));

#ifdef STATS_ON
    /* Collect stats */
    dpu_tasklet_stats_t tasklet_stats[nb_dpus_per_rank][NR_TASKLETS];
    for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
        DPU_FOREACH (rank, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &tasklet_stats[each_dpu][each_tasklet]));
        }
        DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_FROM_DPU, XSTR(DPU_TASKLET_STATS_VAR), each_tasklet * sizeof(dpu_tasklet_stats_t),
            sizeof(dpu_tasklet_stats_t), DPU_XFER_DEFAULT));
    }
#endif

    pthread_mutex_lock(&devices.log_mutex);
    fprintf(devices.log_file, "rank %u offset %u\n", rank_id, dpu_offset);
    unsigned int nb_dpu = get_nb_dpu();
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu >= nb_dpu)
            break;
        fprintf(devices.log_file, "LOG DPU=%u TIME=%llu\n", this_dpu, (unsigned long long)compute_time[each_dpu]);
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
        DPU_ASSERT(dpu_log_read(dpu, devices.log_file));
    }
    fflush(devices.log_file);
    pthread_mutex_unlock(&devices.log_mutex);
}

void dpu_try_get_results_and_log(unsigned int rank_id, unsigned int dpu_offset)
{
    static dpu_result_out_t dummy_results[MAX_DPU_RESULTS];
    static nb_result_t dummy_nb_results;

    struct dpu_set_t rank = devices.ranks[rank_id];
    uint32_t nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    uint8_t *results[nb_dpus_per_rank];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = get_nb_dpu();

    dpu_offset += devices.rank_mram_offset[rank_id];

    dpu_result_out_t *results_buffer[nb_dpus_per_rank];
    nb_result_t *nb_res = get_mem_dpu_nb_res(dpu_offset);
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu >= nb_dpu)
            continue;
        results_buffer[each_dpu] = get_mem_dpu_res(this_dpu);
    }

    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            results[each_dpu] = (uint8_t *)&nb_res[each_dpu];
        } else {
            results[each_dpu] = (uint8_t *)&dummy_nb_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_FROM_DPU, XSTR(DPU_NB_RESULT_VAR), 0, sizeof(nb_result_t), DPU_XFER_DEFAULT));
    unsigned int max_nb_result = 0;
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            results[each_dpu] = (uint8_t *)results_buffer[each_dpu];
            max_nb_result = MAX(max_nb_result, nb_res[each_dpu]);
        } else {
            results[each_dpu] = (uint8_t *)dummy_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(
        rank, DPU_XFER_FROM_DPU, XSTR(DPU_RESULT_VAR), 0, (max_nb_result + 1) * sizeof(dpu_result_out_t), DPU_XFER_DEFAULT));

    uint32_t checksums[nb_dpus_per_rank];
    DPU_FOREACH(rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &checksums[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_FROM_DPU, "output_checksum", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));
    DPU_FOREACH(rank, dpu, each_dpu) {
        uint32_t host_checksum = compute_checksum((uint32_t *)(results[each_dpu]), nb_res[each_dpu] * sizeof(dpu_result_out_t));
        if (checksums[each_dpu] != host_checksum) {
            fprintf(stderr, "ERROR in checksum of dpu %u.%u (DPU#%u) (0x%x != 0x%x): mram_read failed\n", rank_id, each_dpu,
                dpu_offset + each_dpu, checksums[each_dpu], host_checksum);
        }
#ifdef ASSERT_CHECKSUM_REF
        if (checksums[each_dpu]
            != checksums_ref[(dpu_offset + each_dpu) + 640 * (dpus_mgmt_pass[rank_id] + 250 * (dpus_mgmt_round))]) {
            fprintf(stderr, "ERROR in checksum of dpu %u.%u (DPU#%u) (0x%x != 0x%x): compute on dpu failed\n", rank_id, each_dpu,
                dpu_offset + each_dpu, checksums[each_dpu],
                checksums_ref[(dpu_offset + each_dpu) + 640 * (dpus_mgmt_pass[rank_id] + 250 * (dpus_mgmt_round))]);
        }
#endif
#ifdef GENERATE_CHECKSUM_REF
        checksums_ref[(dpu_offset + each_dpu) + 640 * (dpus_mgmt_pass[rank_id] + 250 * (dpus_mgmt_round))] = checksums[each_dpu];
#endif
    }

    dpu_try_log(rank_id, dpu_offset);
}
