/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <semaphore.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <dpu.h>
#include <dpu_log.h>
#include <dpu_management.h>

#include "accumulateread.h"
#include "common.h"
#include "dispatch.h"
#include "dpu_backend.h"
#include "index.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "upvc.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

DPU_INCBIN(upvc_dpu_program, DPU_BINARY);

struct triplet {
    uint32_t rank;
    uint32_t ci;
    uint32_t dpu;
};

typedef struct {
    unsigned int nb_dpus_per_rank[NB_RANKS_MAX];
    unsigned int rank_mram_offset[NB_RANKS_MAX];
    unsigned int nb_dpus;
    struct dpu_set_t all_ranks;
    struct dpu_set_t ranks[NB_RANKS_MAX];
    pthread_mutex_t log_mutex;
    FILE *log_file;
    struct triplet *dpus;
} devices_t;

static bool dpu_backend_initialized = false;
static devices_t devices;

static void dpu_try_write_dispatch_into_mram(unsigned int dpu_offset, unsigned int pass_id)
{
    static const dpu_request_t dummy_dpu_requests[MAX_DPU_REQUEST];
    static dispatch_request_t dummy_dispatch = {
        .nb_reads = 0,
        .dpu_requests = (dpu_request_t *)dummy_dpu_requests,
    };

    dispatch_request_t *io_header[devices.nb_dpus];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = index_get_nb_dpu();
    unsigned int max_dispatch_size = 0;
    DPU_FOREACH (devices.all_ranks, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu < nb_dpu) {
            io_header[each_dpu] = dispatch_get(this_dpu, pass_id);
        } else {
            io_header[each_dpu] = &dummy_dispatch;
        }
        max_dispatch_size = MAX(max_dispatch_size, io_header[each_dpu]->nb_reads * sizeof(dpu_request_t));
    }
    DPU_FOREACH (devices.all_ranks, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &io_header[each_dpu]->nb_reads));
    }
    DPU_ASSERT(dpu_push_xfer(devices.all_ranks, DPU_XFER_TO_DPU, XSTR(DPU_NB_REQUEST_VAR), 0, sizeof(nb_request_t), DPU_XFER_ASYNC));

    DPU_FOREACH (devices.all_ranks, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, io_header[each_dpu]->dpu_requests));
    }
    DPU_ASSERT(dpu_push_xfer(devices.all_ranks, DPU_XFER_TO_DPU, XSTR(DPU_REQUEST_VAR), 0, max_dispatch_size, DPU_XFER_ASYNC));
}

static __attribute__((used)) dpu_error_t dpu_try_log(struct dpu_set_t rank, uint32_t rank_id, void *arg)
{
    unsigned int dpu_offset = (unsigned int)(uintptr_t)arg;
    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    dpu_compute_time_t compute_time[nb_dpus_per_rank];
    struct dpu_set_t dpu;
    unsigned int each_dpu;
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &compute_time[each_dpu]));
    }
    DPU_ASSERT(
        dpu_push_xfer(rank, DPU_XFER_FROM_DPU, XSTR(DPU_COMPUTE_TIME_VAR), 0, sizeof(dpu_compute_time_t), DPU_XFER_DEFAULT));

    /* Collect stats */
    dpu_tasklet_stats_t tasklet_stats[nb_dpus_per_rank][NR_TASKLETS];
    for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
        DPU_FOREACH (rank, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &tasklet_stats[each_dpu][each_tasklet]));
        }
        DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_FROM_DPU, XSTR(DPU_TASKLET_STATS_VAR), each_tasklet * sizeof(dpu_tasklet_stats_t),
            sizeof(dpu_tasklet_stats_t), DPU_XFER_DEFAULT));
    }

    pthread_mutex_lock(&devices.log_mutex);
    fprintf(devices.log_file, "rank %u offset %u\n", rank_id, dpu_offset);
    unsigned int nb_dpu = index_get_nb_dpu();
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu >= nb_dpu)
            break;
        fprintf(devices.log_file, "LOG DPU=%u TIME=%llu\n", this_dpu, (unsigned long long)compute_time[each_dpu]);
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
        fprintf(devices.log_file, "LOG DPU=%u REQ=%u\n", this_dpu, agreagated_stats.nb_reqs);
        fprintf(devices.log_file, "LOG DPU=%u NODP=%u\n", this_dpu, agreagated_stats.nb_nodp_calls);
        fprintf(devices.log_file, "LOG DPU=%u ODPD=%u\n", this_dpu, agreagated_stats.nb_odpd_calls);
        fprintf(devices.log_file, "LOG DPU=%u NODP_TIME=%llu\n", this_dpu, (unsigned long long)agreagated_stats.nodp_time);
        fprintf(devices.log_file, "LOG DPU=%u ODPD_TIME=%llu\n", this_dpu, (unsigned long long)agreagated_stats.odpd_time);
        fprintf(devices.log_file, "LOG DPU=%u RESULTS=%u\n", this_dpu, agreagated_stats.nb_results);
        fprintf(devices.log_file, "LOG DPU=%u DATA_IN=%u\n", this_dpu, agreagated_stats.mram_data_load);
        fprintf(devices.log_file, "LOG DPU=%u RESULT_OUT=%u\n", this_dpu, agreagated_stats.mram_result_store);
        fprintf(devices.log_file, "LOG DPU=%u LOAD=%u\n", this_dpu, agreagated_stats.mram_load);
        fprintf(devices.log_file, "LOG DPU=%u STORE=%u\n", this_dpu, agreagated_stats.mram_store);

        DPU_ASSERT(dpu_log_read(dpu, devices.log_file));
    }
    fflush(devices.log_file);
    pthread_mutex_unlock(&devices.log_mutex);
    return DPU_OK;
}

typedef union {
    struct {
        uint32_t dpu_offset;
        uint32_t pass_id;
    };
    uint64_t info;
} pass_info_t;

_Static_assert(sizeof(pass_info_t) == sizeof(uint64_t), "dpu_callback using this type will not be functional");

static dpu_error_t dpu_get_results(struct dpu_set_t rank, uint32_t rank_id, __attribute__((unused)) void *arg)
{
    static dpu_result_out_t dummy_results[MAX_DPU_RESULTS];
    pass_info_t info = (pass_info_t)(uintptr_t)arg;
    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = index_get_nb_dpu();
    unsigned int max_nb_result = 0;
    unsigned int mram_offset = devices.rank_mram_offset[rank_id];
    unsigned int dpu_offset = info.dpu_offset + mram_offset;
    unsigned int pass_id = info.pass_id;

    DPU_FOREACH (rank, dpu, each_dpu) {
        uint8_t *results;
        if ((each_dpu + dpu_offset) < nb_dpu) {
            acc_results_t *acc_res = accumulate_get_buffer(each_dpu + mram_offset, pass_id);
            results = (uint8_t *)acc_res->results;
            max_nb_result = MAX(max_nb_result, acc_res->nb_res);
        } else {
            results = (uint8_t *)dummy_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results));
    }
    DPU_ASSERT(dpu_push_xfer(
        rank, DPU_XFER_FROM_DPU, XSTR(DPU_RESULT_VAR), 0, (max_nb_result + 1) * sizeof(dpu_result_out_t), DPU_XFER_DEFAULT));

    return DPU_OK;
}

static void dpu_try_get_results_and_log(unsigned int dpu_offset, unsigned int pass_id)
{
    static nb_result_t dummy_nb_results;

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = index_get_nb_dpu();
    DPU_FOREACH (devices.all_ranks, dpu, each_dpu) {
        uint8_t *results;
        if ((each_dpu + dpu_offset) < nb_dpu) {
            results = (uint8_t *)&(accumulate_get_buffer(each_dpu, pass_id)->nb_res);
        } else {
            results = (uint8_t *)&dummy_nb_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results));
    }
    DPU_ASSERT(dpu_push_xfer(devices.all_ranks, DPU_XFER_FROM_DPU, XSTR(DPU_NB_RESULT_VAR), 0, sizeof(nb_result_t), DPU_XFER_ASYNC));

    pass_info_t info = { .dpu_offset = dpu_offset, .pass_id = pass_id };
    DPU_ASSERT(dpu_callback(devices.all_ranks, dpu_get_results, (void *)info.info, DPU_CALLBACK_ASYNC));
#ifdef STATS_ON
    DPU_ASSERT(dpu_callback(devices.all_ranks, dpu_try_log, (void *)dpu_offset, DPU_CALLBACK_ASYNC));
#endif
}

static dpu_error_t sem_post_dispatch_free_sem(
    __attribute__((unused)) struct dpu_set_t set, __attribute__((unused)) uint32_t rank_id, void *arg)
{
    sem_post((sem_t *)arg);
    return DPU_OK;
}

static dpu_error_t sem_post_exec_to_acc_sem(
    __attribute__((unused)) struct dpu_set_t set, __attribute__((unused)) uint32_t rank_id, void *arg)
{
    sem_post((sem_t *)arg);
    return DPU_OK;
}

void run_on_dpu(unsigned int dpu_offset, unsigned int pass_id, sem_t *dispatch_free_sem, sem_t *acc_wait_sem,
    sem_t *exec_to_acc_sem, sem_t *dispatch_to_exec_sem)
{
    dpu_try_write_dispatch_into_mram(dpu_offset, pass_id);

    DPU_ASSERT(dpu_callback(devices.all_ranks, sem_post_dispatch_free_sem, dispatch_free_sem,
        DPU_CALLBACK_ASYNC | DPU_CALLBACK_NONBLOCKING | DPU_CALLBACK_SINGLE_CALL));
    DPU_ASSERT(dpu_launch(devices.all_ranks, DPU_ASYNCHRONOUS));
    sem_wait(acc_wait_sem);

    dpu_try_get_results_and_log(dpu_offset, pass_id);

    DPU_ASSERT(dpu_callback(devices.all_ranks, sem_post_exec_to_acc_sem, exec_to_acc_sem,
        DPU_CALLBACK_ASYNC | DPU_CALLBACK_NONBLOCKING | DPU_CALLBACK_SINGLE_CALL));
    sem_wait(dispatch_to_exec_sem);
}

void init_backend_dpu(unsigned int *nb_dpus_per_run)
{
    const char *profile = "cycleAccurate=true,nrJobsPerRank=64";

    DPU_ASSERT(dpu_alloc(get_nb_dpu(), profile, &devices.all_ranks));
    DPU_ASSERT(dpu_load_from_incbin(devices.all_ranks, &upvc_dpu_program, NULL));
    DPU_ASSERT(dpu_get_nr_dpus(devices.all_ranks, nb_dpus_per_run));

    unsigned int each_rank;
    struct dpu_set_t rank;
    devices.nb_dpus = 0;
    DPU_RANK_FOREACH (devices.all_ranks, rank, each_rank) {
        devices.ranks[each_rank] = rank;
        DPU_ASSERT(dpu_get_nr_dpus(rank, &devices.nb_dpus_per_rank[each_rank]));
        devices.rank_mram_offset[each_rank] = devices.nb_dpus;
        devices.nb_dpus += devices.nb_dpus_per_rank[each_rank];
    }
    printf("%u DPUs allocated\n", devices.nb_dpus);
    assert(devices.nb_dpus == *nb_dpus_per_run);

#ifdef STATS_ON
    pthread_mutex_init(&devices.log_mutex, NULL);
    char filename[1024];
    sprintf(filename, "%s_log.txt", get_input_path());
    devices.log_file = fopen(filename, "w");
    CHECK_FILE(devices.log_file, filename);
#endif

    struct dpu_set_t dpu;
    uint32_t each_dpu;
    devices.dpus = malloc(sizeof(struct triplet) * devices.nb_dpus);
    DPU_FOREACH (devices.all_ranks, dpu, each_dpu) {
        struct dpu_t *dpu_t = dpu_from_set(dpu);
        devices.dpus[each_dpu].rank = dpu_get_rank_id(dpu_get_rank(dpu_t));
        devices.dpus[each_dpu].ci = dpu_get_slice_id(dpu_t);
        devices.dpus[each_dpu].dpu = dpu_get_member_id(dpu_t);
    }
    dpu_backend_initialized = true;
}

void free_backend_dpu()
{
    DPU_ASSERT(dpu_free(devices.all_ranks));
#ifdef STATS_ON
    pthread_mutex_destroy(&devices.log_mutex);
    fclose(devices.log_file);
#endif
}

typedef union {
    struct {
        uint32_t dpu_offset;
        uint32_t delta_neighbour;
    };
    uint64_t info;
} load_info_t;

static dpu_error_t load_mram_rank(struct dpu_set_t rank, uint32_t rank_id, void *args) {
    load_info_t info = (load_info_t)(uintptr_t)args;
    unsigned int dpu_offset = info.dpu_offset;
    unsigned int delta_neighbour = info.delta_neighbour;
    unsigned int nb_dpu = index_get_nb_dpu();
    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    uint8_t *mram[nb_dpus_per_rank];
    size_t mram_size[nb_dpus_per_rank];

    memset(mram, 0, sizeof(uint8_t *) * nb_dpus_per_rank);
    memset(mram_size, 0, sizeof(size_t) * nb_dpus_per_rank);

    unsigned int each_dpu;
    struct dpu_set_t dpu;
    dpu_offset += devices.rank_mram_offset[rank_id];
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = dpu_offset + each_dpu;
        if (this_dpu < nb_dpu) {
            mram_size[each_dpu] = mram_load(&mram[each_dpu], this_dpu);
        }
    }

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
    return DPU_OK;
}

void load_mram_dpu(unsigned int dpu_offset, int delta_neighbour)
{
    load_info_t info = {.dpu_offset = dpu_offset, delta_neighbour = delta_neighbour};
    dpu_callback(devices.all_ranks, load_mram_rank, (void *)info.info, DPU_CALLBACK_ASYNC);
}

void wait_dpu_dpu() {
    DPU_ASSERT(dpu_sync(devices.all_ranks));
}

void get_dpu_info(uint32_t numdpu, uint32_t *rank, uint32_t *ci, uint32_t *dpu)
{
    if (!dpu_backend_initialized) {
        return;
    }
    *rank = devices.dpus[numdpu].rank;
    *ci = devices.dpus[numdpu].ci;
    *dpu = devices.dpus[numdpu].dpu;
}
