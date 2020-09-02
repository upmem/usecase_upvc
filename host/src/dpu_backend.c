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
#include <dpu_runner.h>

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

typedef struct {
    unsigned int nb_dpus_per_rank[NB_RANKS_MAX];
    unsigned int rank_mram_offset[NB_RANKS_MAX];
    struct dpu_set_t all_ranks;
    struct dpu_set_t ranks[NB_RANKS_MAX];
    pthread_mutex_t log_mutex;
    FILE *log_file;
} devices_t;

static devices_t devices;

void get_dpu_id(unsigned int numdpu, unsigned int *rank_id, unsigned *slice_id, unsigned *dpu_id) {
  unsigned int i = 0;
  struct dpu_set_t dpu;
  uint32_t each_dpu;
  uint32_t nb_ranks;
  dpu_get_nr_ranks(devices.all_ranks, &nb_ranks);
  while (i + 1 < nb_ranks && numdpu > devices.rank_mram_offset[i + 1]) i++;

  *rank_id = -1;
  *slice_id = -1;
  *dpu_id = -1;
  DPU_FOREACH(devices.ranks[i], dpu, each_dpu) {
    if (each_dpu == numdpu - devices.rank_mram_offset[i]) {
      *rank_id = dpu_get_rank_id(dpu_get_rank(dpu_from_set(dpu)));
      *slice_id = dpu_get_slice_id(dpu_from_set(dpu));
      *dpu_id = dpu_get_member_id(dpu_from_set(dpu));
    }
  }
}

static void dpu_try_run(unsigned int rank_id)
{
    struct dpu_set_t rank = devices.ranks[rank_id];
    uint32_t nb_dpu_running;
    int timeout = 5000000;
    dpu_error_t status = dpu_boot_rank(dpu_rank_from_set(rank)); //dpu_launch(rank, DPU_ASYNCHRONOUS);
    uint32_t ret;

    do {
        ret = dpu_poll_rank(dpu_rank_from_set(rank));
        if (ret != DPU_OK){
            printf("ERROR polling rank %d\n", rank_id);
            return;
        }
        nb_dpu_running = dpu_get_run_context(dpu_rank_from_set(rank))->nb_dpu_running;
    } while (nb_dpu_running && timeout-- > 0);
    if (timeout <= 0) {
        printf("ERROR timeout polling rank %d\n", rank_id);
        return;
    }

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

const char *get_nb_result_var(unsigned int pass_id) {
  unsigned int pass_id_mod = pass_id % 4;
    if (pass_id_mod == 0)
        return "DPU_NB_RESULT_VAR0";
    else if (pass_id_mod == 1)
        return "DPU_NB_RESULT_VAR1";
    else if (pass_id_mod == 2)
      return "DPU_NB_RESULT_VAR2";
    else
      return "DPU_NB_RESULT_VAR3";
}
const char *get_result_var(unsigned int pass_id) {
  unsigned int pass_id_mod = pass_id % 4;
    if (pass_id_mod == 0)
        return "DPU_RESULT_VAR0";
    else if (pass_id_mod == 1)
        return "DPU_RESULT_VAR1";
    else if (pass_id_mod == 2)
      return "DPU_RESULT_VAR2";
    else
      return "DPU_RESULT_VAR3";
}
const char *get_nb_request_var(unsigned int pass_id) {
  unsigned int pass_id_mod = pass_id % 4;
    if (pass_id_mod == 0)
        return "DPU_NB_REQUEST_VAR0";
    else if (pass_id_mod == 1)
        return "DPU_NB_REQUEST_VAR1";
    else if (pass_id_mod == 2)
      return "DPU_NB_REQUEST_VAR2";
    else
      return "DPU_NB_REQUEST_VAR3";
}
const char *get_request_var(unsigned int pass_id) {
  unsigned int pass_id_mod = pass_id % 4;
    if (pass_id_mod == 0)
        return "DPU_REQUEST_VAR0";
    else if (pass_id_mod == 1)
        return "DPU_REQUEST_VAR1";
    else if (pass_id_mod == 2)
      return "DPU_REQUEST_VAR2";
    else
      return "DPU_REQUEST_VAR3";
}

static void dpu_try_write_dispatch_into_mram(unsigned int rank_id, unsigned int dpu_offset, unsigned int pass_id)
{
    static const dpu_request_t dummy_dpu_requests[MAX_DPU_REQUEST];
    static dispatch_request_t dummy_dispatch = {
        .nb_reads = 0,
        .dpu_requests = (dpu_request_t *)dummy_dpu_requests,
    };
    struct dpu_set_t rank = devices.ranks[rank_id];

    unsigned int nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    dispatch_request_t *io_header[nb_dpus_per_rank];

    uint32_t dpu_result_var_idx = pass_id % 4;
    DPU_ASSERT(dpu_copy_to(rank, "dpu_result_var_idx", 0, &dpu_result_var_idx, sizeof(uint32_t)));

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = index_get_nb_dpu();
    unsigned int max_dispatch_size = 0;
    dpu_offset += devices.rank_mram_offset[rank_id];
    DPU_FOREACH (rank, dpu, each_dpu) {
        unsigned int this_dpu = each_dpu + dpu_offset;
        if (this_dpu < nb_dpu) {
            io_header[each_dpu] = dispatch_get(this_dpu, pass_id);
        } else {
            io_header[each_dpu] = &dummy_dispatch;
        }
        max_dispatch_size = MAX(max_dispatch_size, io_header[each_dpu]->nb_reads * sizeof(dpu_request_t));
    }
    const char *dpu_nb_request_var = get_nb_request_var(pass_id);
    const char *dpu_request_var = get_request_var(pass_id);
    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &io_header[each_dpu]->nb_reads));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, dpu_nb_request_var, 0, sizeof(nb_request_t), DPU_XFER_DEFAULT));

    DPU_FOREACH (rank, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, io_header[each_dpu]->dpu_requests));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_TO_DPU, dpu_request_var, 0, max_dispatch_size, DPU_XFER_DEFAULT));
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
    unsigned int nb_dpu = index_get_nb_dpu();
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

#endif
        DPU_ASSERT(dpu_log_read(dpu, devices.log_file));
    }
    fflush(devices.log_file);
    pthread_mutex_unlock(&devices.log_mutex);
}

static void dpu_try_get_results_and_log(unsigned int rank_id, unsigned int dpu_offset, unsigned int pass_id)
{
    static dpu_result_out_t dummy_results[MAX_DPU_RESULTS];
    static nb_result_t dummy_nb_results;

    struct dpu_set_t rank = devices.ranks[rank_id];
    uint32_t nb_dpus_per_rank = devices.nb_dpus_per_rank[rank_id];
    uint8_t *results[nb_dpus_per_rank];

    struct dpu_set_t dpu;
    unsigned int each_dpu;
    unsigned int nb_dpu = index_get_nb_dpu();

    unsigned int mram_offset = devices.rank_mram_offset[rank_id];
    dpu_offset += mram_offset;

    const char *dpu_nb_result_var = get_nb_result_var(pass_id);
    const char *dpu_result_var = get_result_var(pass_id);

    DPU_FOREACH (rank, dpu, each_dpu) {
        if ((each_dpu + dpu_offset) < nb_dpu) {
            results[each_dpu] = (uint8_t *)&(accumulate_get_buffer(each_dpu + mram_offset, pass_id)->nb_res);
        } else {
            results[each_dpu] = (uint8_t *)&dummy_nb_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(rank, DPU_XFER_FROM_DPU, dpu_nb_result_var, 0, sizeof(nb_result_t), DPU_XFER_DEFAULT));

    unsigned int max_nb_result = 0;
    DPU_FOREACH (rank, dpu, each_dpu) {
        if ((each_dpu + dpu_offset) < nb_dpu) {
            acc_results_t *acc_res = accumulate_get_buffer(each_dpu + mram_offset, pass_id);
            results[each_dpu] = (uint8_t *)acc_res->results;
            max_nb_result = MAX(max_nb_result, acc_res->nb_res);
        } else {
            results[each_dpu] = (uint8_t *)dummy_results;
        }
        DPU_ASSERT(dpu_prepare_xfer(dpu, results[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(
        rank, DPU_XFER_FROM_DPU, dpu_result_var, 0, (max_nb_result + 1) * sizeof(dpu_result_out_t), DPU_XFER_DEFAULT));

    dpu_try_log(rank_id, dpu_offset);
}

void run_on_dpu(
    unsigned int dpu_offset, unsigned int rank_id, unsigned int pass_id, sem_t *dispatch_free_sem, sem_t *acc_wait_sem)
{
    double t1, t2, t3, t4, t5;
    PRINT_TIME_WRITE_READS(rank_id);
    t1 = my_clock();
    dpu_try_write_dispatch_into_mram(rank_id, dpu_offset, pass_id);
    sem_post(dispatch_free_sem);

    t2 = my_clock();
    PRINT_TIME_WRITE_READS(rank_id);
    PRINT_TIME_COMPUTE(rank_id);

    dpu_try_run(rank_id);

    PRINT_TIME_COMPUTE(rank_id);
    t3 = my_clock();

    sem_wait(acc_wait_sem);

    t4 = my_clock();
    PRINT_TIME_READ_RES(rank_id);
    dpu_try_get_results_and_log(rank_id, dpu_offset, pass_id);
    PRINT_TIME_READ_RES(rank_id);
    t5 = my_clock();

    print(rank_id, "W %.3lf", t2 - t1);
    print(rank_id, "C %.3lf", t3 - t2);
    print(rank_id, "R %.3lf", t5 - t4);
}

void init_backend_dpu(unsigned int *nb_dpus_per_run, unsigned int *nb_ranks_per_run)
{
    const char *profile = "cycleAccurate=true,nrThreadsPerRank=0";

    DPU_ASSERT(dpu_alloc(get_nb_dpu(), profile, &devices.all_ranks));
    DPU_ASSERT(dpu_load_from_incbin(devices.all_ranks, &upvc_dpu_program, NULL));
    DPU_ASSERT(dpu_get_nr_ranks(devices.all_ranks, nb_ranks_per_run));
    DPU_ASSERT(dpu_get_nr_dpus(devices.all_ranks, nb_dpus_per_run));
    assert(*nb_ranks_per_run <= NB_RANKS_MAX);

    unsigned int nb_dpus = 0;
    unsigned int each_rank;
    struct dpu_set_t rank;
    DPU_RANK_FOREACH (devices.all_ranks, rank, each_rank) {
        devices.ranks[each_rank] = rank;
        DPU_ASSERT(dpu_get_nr_dpus(rank, &devices.nb_dpus_per_rank[each_rank]));
        devices.rank_mram_offset[each_rank] = nb_dpus;
        nb_dpus += devices.nb_dpus_per_rank[each_rank];
    }
    assert(nb_dpus == *nb_dpus_per_run);
    printf("%u DPUs allocated\n", nb_dpus);
    unsigned int rank_id, slice_id, dpu_id;
    get_dpu_id(1206, &rank_id, &slice_id, &dpu_id);
    printf("%d %u.%u.%u\n", 356, rank_id, slice_id, dpu_id);
    pthread_mutex_init(&devices.log_mutex, NULL);
    char filename[1024];
    sprintf(filename, "%s_log.txt", get_input_path());
    devices.log_file = fopen(filename, "w");
    CHECK_FILE(devices.log_file, filename);
}

void free_backend_dpu()
{
    DPU_ASSERT(dpu_free(devices.all_ranks));
    pthread_mutex_destroy(&devices.log_mutex);
    fclose(devices.log_file);
}

void load_mram_dpu(unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour)
{
    unsigned int nb_dpu = index_get_nb_dpu();
    struct dpu_set_t rank = devices.ranks[rank_id];
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
}

int get_rank_numa_node(int rank_id) {
    return dpu_get_numa_node(dpu_rank_from_set(devices.ranks[rank_id]));
}

