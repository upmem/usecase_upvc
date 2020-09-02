/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdint.h>
#include <string.h>

#include <alloc.h>
#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <perfcounter.h>

#include "debug.h"
#include "dout.h"
#include "nodp.h"
#include "odpd.h"
#include "request_pool.h"
#include "result_pool.h"
#include "stats.h"

#include "common.h"

BARRIER_INIT(init_barrier, NR_TASKLETS);

STDOUT_BUFFER_INIT(128 * 1024);

/**
 * @brief The MRAM information, shared between the tasklets.
 *
 * Information are gathered by the first tasklet during boot and used by other tasklets
 * to get the memory topology.
 */
__host delta_info_t DPU_MRAM_INFO_VAR;

/**
 * @brief The statistic of each tasklet in the mram.
 */
__mram_noinit dpu_tasklet_stats_t DPU_TASKLET_STATS_VAR[NR_TASKLETS];
#ifdef STATS_ON
#define DPU_TASKLET_STATS_WRITE(res, addr)                                                                                       \
    do {                                                                                                                         \
        mram_write(res, addr, sizeof(dpu_tasklet_stats_t));                                                                      \
    } while (0)
#else
#define DPU_TASKLET_STATS_WRITE(res, addr)
#endif

/**
 * @brief The compute time performance value store in mram.
 */
__host dpu_compute_time_t DPU_COMPUTE_TIME_VAR;

/**
 * @brief Maximum score allowed.
 */
#define MAX_SCORE (40)

/**
 * @brief Number of reference read to be fetch per mram read
 */
#define NB_REF_PER_READ (4)

/**
 * @brief Global table of dout_t structure.
 *
 * As the structure is quite big, let's have it in the heap and not in the stack.
 * Each tasklet will use the dout_t at the index of its tasklet_id.
 * It uses this structure to store temporary local results.
 */
__dma_aligned static dout_t global_dout[NR_TASKLETS];

__dma_aligned coords_and_nbr_t coords_and_nbr[NR_TASKLETS][NB_REF_PER_READ];
__dma_aligned dpu_request_t requests[NR_TASKLETS];

/**
 * @brief Fetches a neighbour from the neighbour area.
 *
 * @param base    Offset to the first neighbour of the requested pool within the neighbour area.
 * @param idx     Index of this neighbour in the specified pool.
 * @param buffer  To contain the result, must be the size of a neighbour and coordinates, plus 64 bits alignment.
 * @param stats   To update statistical report.
 */
static void load_reference_multiple_nbr_and_coords_at(
    unsigned int base, unsigned int idx, uint8_t *cache, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    /* The input starts with coordinates (8 bytes), followed by the neighbour. Structure is aligned
     * on 8 bytes boundary.
     */
    uintptr_t coords_nbr_address = ((uintptr_t)DPU_MRAM_HEAP_POINTER + (base + idx) * sizeof(coords_and_nbr_t));
    unsigned int coords_nbr_len_total = sizeof(coords_and_nbr_t) * NB_REF_PER_READ;
    ASSERT_DMA_ADDR(coords_nbr_address, cache, coords_nbr_len_total);
    ASSERT_DMA_LEN(coords_nbr_len_total);
    mram_read((__mram_ptr void *)coords_nbr_address, cache, coords_nbr_len_total);
    STATS_INCR_LOAD(stats, coords_nbr_len_total);
    STATS_INCR_LOAD_DATA(stats, coords_nbr_len_total);
}

static void get_time_and_accumulate(dpu_compute_time_t *accumulate_time, perfcounter_t *last_time)
{
    perfcounter_t current_time = perfcounter_get();
    if (current_time < *last_time) {
        *accumulate_time = *accumulate_time + 0x1000000000;
    }
    *last_time = current_time;
}

static void compare_neighbours(sysname_t tasklet_id, uint32_t *mini, coords_and_nbr_t *cached_coords_and_nbr,
    uint8_t *current_read_nbr, dpu_request_t *request, dout_t *dout, dpu_tasklet_stats_t *tasklet_stats)
{
    uint32_t score, score_nodp, score_odpd = UINT_MAX;
    uint8_t *ref_nbr = &cached_coords_and_nbr->nbr[0];
    STATS_TIME_VAR(start, end, acc);

    STATS_GET_START_TIME(start, acc, end);

    score = score_nodp = nodp(current_read_nbr, ref_nbr, *mini, SIZE_NEIGHBOUR_IN_BYTES - DPU_MRAM_INFO_VAR);

    STATS_GET_END_TIME(end, acc);
    STATS_STORE_NODP_TIME(tasklet_stats, (end + acc - start));
    STATS_INCR_NB_NODP_CALLS(*tasklet_stats);

    if (score_nodp == UINT_MAX) {
        STATS_GET_START_TIME(start, acc, end);

        score_odpd = score = odpd(current_read_nbr, ref_nbr, *mini, NB_BYTES_TO_SYMS(SIZE_NEIGHBOUR_IN_BYTES, DPU_MRAM_INFO_VAR));

        STATS_GET_END_TIME(end, acc);
        STATS_STORE_ODPD_TIME(tasklet_stats, (end + acc - start));
        STATS_INCR_NB_ODPD_CALLS(*tasklet_stats);
    }

    if (score > *mini) {
        return;
    }

    if (score < *mini) {
        *mini = score;
        /* Get rid of previous results, we found a better one */
        dout_clear(dout);
    }

    if (dout->nb_results >= MAX_RESULTS_PER_READ) {
        printf("WARNING! too many results for request!\n");
        /* Trigger a fault, since this should never happen. */
        halt();
    }

    dout_add(dout, request->num, (unsigned int)score, cached_coords_and_nbr->coord.seed_nr, cached_coords_and_nbr->coord.seq_nr,
        tasklet_stats);
}

static void compute_request(sysname_t tasklet_id, coords_and_nbr_t *cached_coords_and_nbr, uint8_t *current_read_nbr,
    dpu_request_t *request, dout_t *dout, dpu_tasklet_stats_t *tasklet_stats)
{
    uint32_t mini = MAX_SCORE;
    for (unsigned int idx = 0; idx < request->count; idx += NB_REF_PER_READ) {
        load_reference_multiple_nbr_and_coords_at(request->offset, idx, (uint8_t *)cached_coords_and_nbr, tasklet_stats);
        for (unsigned int ref_id = 0; ref_id < NB_REF_PER_READ && ((idx + ref_id) < request->count); ref_id++) {
            compare_neighbours(tasklet_id, &mini, &cached_coords_and_nbr[ref_id], current_read_nbr, request, dout, tasklet_stats);
        }
    }
}

/**
 * @brief Executes the mapping/align procedure.
 *
 * This function is executed by every tasklet, which picks the requests from the request pool
 * and perform the comparison with the reference genome.
 *
 * @return zero.
 */
static void run_align(sysname_t tasklet_id, dpu_compute_time_t *accumulate_time, perfcounter_t *current_time)
{
    STATS_ATTRIBUTE __dma_aligned dpu_tasklet_stats_t tasklet_stats = {
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
    dout_t *dout = &global_dout[tasklet_id];
    coords_and_nbr_t *cached_coords_and_nbr = coords_and_nbr[tasklet_id];
    dpu_request_t *request = &requests[tasklet_id];
    uint8_t *current_read_nbr = &request->nbr[0];

    dout_init(tasklet_id, dout);

    while (request_pool_next(request, &tasklet_stats)) {
        if (tasklet_id == 0) {
            get_time_and_accumulate(accumulate_time, current_time);
        }

        STATS_INCR_NB_REQS(tasklet_stats);

        dout_clear(dout);

        compute_request(tasklet_id, cached_coords_and_nbr, current_read_nbr, request, dout, &tasklet_stats);

        STATS_INCR_NB_RESULTS(tasklet_stats, dout->nb_results);
        result_pool_write(dout, &tasklet_stats);
    }

    DPU_TASKLET_STATS_WRITE(&tasklet_stats, (__mram_ptr void *)(&DPU_TASKLET_STATS_VAR[tasklet_id]));

    result_pool_finish(&tasklet_stats);
}

/**
 * @brief Common main: gather the seed mapping for this DPU, then start processing the requests on every tasklet.
 *
 * @return 0 if everything worked fine
 */
int main()
{
    sysname_t tasklet_id = me();
    perfcounter_t start_time, current_time;

    if (tasklet_id == 0) {
        DPU_COMPUTE_TIME_VAR = 0ULL;
        mem_reset();
        perfcounter_config(COUNT_CYCLES, true);
        current_time = start_time = perfcounter_get();

        check_end_mark();

        request_pool_init();
        result_pool_init();

        if (((NB_BYTES_TO_SYMS(SIZE_NEIGHBOUR_IN_BYTES, DPU_MRAM_INFO_VAR) + 2) * 3 * 16) >= 0x10000) {
            printf("cannot run code: symbol length is larger than mulub operation\n");
            halt();
        }
    }

    barrier_wait(&init_barrier);

    run_align(tasklet_id, &DPU_COMPUTE_TIME_VAR, &current_time);

    barrier_wait(&init_barrier);

    if (tasklet_id == 0) {
        check_end_mark();

        get_time_and_accumulate(&DPU_COMPUTE_TIME_VAR, &current_time);
        DPU_COMPUTE_TIME_VAR = (DPU_COMPUTE_TIME_VAR + current_time) - start_time;
    }

    return 0;
}
