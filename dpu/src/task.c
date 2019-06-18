/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <stdint.h>
#include <string.h>

#include <alloc.h>
#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <perfcounter.h>

#define NB_BYTES_TO_SYMS(len, delta) (((len) - (delta)) << 2)

#include "debug.h"
#include "dout.h"
#include "nodp.h"
#include "odpd.h"
#include "request_pool.h"
#include "result_pool.h"
#include "stats.h"

#include "common.h"

MUTEX_INIT(miscellaneous_mutex);
BARRIER_INIT(init_barrier, NB_RUNNING_TASKLETS);

#define TASKLETS_INITIALIZER TASKLETS(16, main, 256, 0)
_Static_assert(16 == NB_TASKLET_PER_DPU, "NB_TASKLET_PER_DPU does not match the number of tasklets declare in rt.h");
#include <rt.h>

/**
 * @brief The statistic of each tasklet in the mram.
 */
__mram dpu_tasklet_stats_t DPU_TASKLET_STATS_VAR[NB_TASKLET_PER_DPU];
#ifdef STATS_ON
#define DPU_TASKLET_STATS_WRITE(res, addr)                                                                                       \
    do {                                                                                                                         \
        mram_write48(res, addr);                                                                                                 \
    } while (0)
#else
#define DPU_TASKLET_STATS_WRITE(res, addr)
#endif
_Static_assert(sizeof(dpu_tasklet_stats_t) == 48,
    "dpu_tasklet_stats_t size changed (make sure that DPU_TASKLET_STATS_WRITE changed as well)");

/**
 * @brief The compute time performance value store in mram.
 */
__mram dpu_compute_time_t DPU_COMPUTE_TIME_VAR;
#define DPU_COMPUTE_TIME_WRITE(res)                                                                                              \
    do {                                                                                                                         \
        mram_write8(res, (mram_addr_t)&DPU_COMPUTE_TIME_VAR);                                                                    \
    } while (0)
_Static_assert(
    sizeof(dpu_compute_time_t) == 8, "dpu_compute_time_t size changed (make sure that DPU_COMPUTE_TIME_WRITE changed as well)");

/**
 * @brief Maximum score allowed.
 */
#define MAX_SCORE (40)

/**
 * @brief Number of reference read to be fetch per mram read
 */
#define NB_REF_PER_READ (8)

/**
 * @brief The MRAM information, shared between the tasklets.
 *
 * Information are gathered by the first tasklet during boot and used by other tasklets
 * to get the memory topology.
 */
__attribute__((aligned(8))) static mram_info_t mram_info;

/**
 * @brief Global table of dout_t structure.
 *
 * As the structure is quite big, let's have it in the heap and not in the stack.
 * Each tasklet will use the dout_t at the index of its tasklet_id.
 * It uses this structure to store temporary local results.
 */
__attribute__((aligned(8))) static dout_t global_dout[NB_TASKLET_PER_DPU];

/**
 * @brief Fetches a neighbour from the neighbour area.
 *
 * @param base    Offset to the first neighbour of the requested pool within the neighbour area.
 * @param idx     Index of this neighbour in the specified pool.
 * @param buffer  To contain the result, must be the size of a neighbour and coordinates, plus 64 bits alignment.
 * @param stats   To update statistical report.
 */
static void load_reference_multiple_nbr_and_coords_at(
    unsigned int base, unsigned int idx, unsigned int coords_nbr_len, uint8_t *cache, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    /* The input starts with coordinates (8 bytes), followed by the neighbour. Structure is aligned
     * on 8 bytes boundary.
     */
    mram_addr_t coords_nbr_address = (mram_addr_t)(DPU_INPUTS_ADDR(DPU_MRAM_HEAP_POINTER) + (base + idx) * coords_nbr_len);
    unsigned int coords_nbr_len_total = coords_nbr_len * NB_REF_PER_READ;
    ASSERT_DMA_ADDR(coords_nbr_address, cache, coords_nbr_len_total);
    ASSERT_DMA_LEN(coords_nbr_len_total);
    mram_readX(coords_nbr_address, cache, coords_nbr_len_total);
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

static void compare_neighbours(sysname_t tasklet_id, int *mini, uint8_t *cached_coords_and_nbr, uint8_t *current_read_nbr,
    dpu_request_t *request, dout_t *dout, dpu_tasklet_stats_t *tasklet_stats, mutex_id_t *mutex_miscellaneous)
{
    int score, score_nodp, score_odpd = -1;
    uint8_t *ref_nbr = cached_coords_and_nbr + sizeof(dpu_result_coord_t);
    STATS_TIME_VAR(start, end, acc);
    DEBUG_RESULTS_VAR;
    DEBUG_REQUESTS_VAR;
    DEBUG_REQUESTS_PRINT_REF(ref_nbr, mram_info.nbr_len, mram_info.delta);

    STATS_GET_START_TIME(start, acc, end);

    score = score_nodp = noDP(current_read_nbr, ref_nbr, mram_info.nbr_len, mram_info.delta, *mini);

    STATS_GET_END_TIME(end, acc);
    STATS_STORE_NODP_TIME(tasklet_stats, (end + acc - start));
    STATS_INCR_NB_NODP_CALLS(*tasklet_stats);

    if (score_nodp == -1) {
        STATS_GET_START_TIME(start, acc, end);

        score_odpd = score
            = odpd(current_read_nbr, ref_nbr, *mini, NB_BYTES_TO_SYMS(mram_info.nbr_len, mram_info.delta), tasklet_id);

        STATS_GET_END_TIME(end, acc);
        STATS_STORE_ODPD_TIME(tasklet_stats, (end + acc - start));
        STATS_INCR_NB_ODPD_CALLS(*tasklet_stats);
    }
    DEBUG_PROCESS_SCORES(mram_info.nbr_len, current_read_nbr, ref_nbr, *mini, score_nodp, score_odpd, *mutex_miscellaneous);

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

    dout_add(dout, request->num, (unsigned int)score, ((uint32_t *)cached_coords_and_nbr)[0],
        ((uint32_t *)cached_coords_and_nbr)[1], tasklet_stats);
    DEBUG_RESULTS_PRINT(tasklet_id, request->num, request->offset, ((uint32_t *)cached_coords_and_nbr)[0],
        ((uint32_t *)cached_coords_and_nbr)[1], score);
}

static void compute_request(sysname_t tasklet_id, unsigned int coords_nbr_len, uint8_t *cached_coords_and_nbr,
    uint8_t *current_read_nbr, dpu_request_t *request, dout_t *dout, dpu_tasklet_stats_t *tasklet_stats,
    mutex_id_t *mutex_miscellaneous)
{
    int mini = MAX_SCORE;
    for (unsigned int idx = 0; idx < request->count; idx += NB_REF_PER_READ) {
        load_reference_multiple_nbr_and_coords_at(request->offset, idx, coords_nbr_len, cached_coords_and_nbr, tasklet_stats);
        for (unsigned int ref_id = 0; ref_id < NB_REF_PER_READ && ((idx + ref_id) < request->count); ref_id++) {
            compare_neighbours(tasklet_id, &mini, cached_coords_and_nbr + (ref_id * coords_nbr_len), current_read_nbr, request,
                dout, tasklet_stats, mutex_miscellaneous);
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
    STATS_ATTRIBUTE __attribute__((aligned(8))) dpu_tasklet_stats_t tasklet_stats = {
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
    mutex_id_t mutex_miscellaneous = MUTEX_GET(miscellaneous_mutex);
    dout_t *dout = &global_dout[tasklet_id];
    unsigned int nbr_len_aligned = ALIGN_DPU(mram_info.nbr_len);
    unsigned int coords_nbr_len = sizeof(dpu_result_coord_t) + nbr_len_aligned;
    uint8_t *cached_coords_and_nbr = mem_alloc_dma(coords_nbr_len * NB_REF_PER_READ);
    uint8_t *request_buffer = mem_alloc_dma(sizeof(dpu_request_t) + nbr_len_aligned);
    uint8_t *current_read_nbr = request_buffer + sizeof(dpu_request_t);
    dpu_request_t *request = (dpu_request_t *)request_buffer;

    DEBUG_RESULTS_VAR;
    DEBUG_REQUESTS_VAR;

    DEBUG_PROCESS_PROFILE(mram_info.nbr_len);

    dout_init(tasklet_id, dout);

    while (request_pool_next(request_buffer, &tasklet_stats)) {
        if (tasklet_id == 0) {
            get_time_and_accumulate(accumulate_time, current_time);
        }

        DEBUG_REQUESTS_PRINT(request, mram_info.nbr_len, mram_info.delta);
        DEBUG_RESULTS_INCR;

        STATS_INCR_NB_REQS(tasklet_stats);

        dout_clear(dout);

        compute_request(tasklet_id, coords_nbr_len, cached_coords_and_nbr, current_read_nbr, request, dout, &tasklet_stats,
            &mutex_miscellaneous);

        STATS_INCR_NB_RESULTS(tasklet_stats, dout->nb_results);
        result_pool_write(dout, &tasklet_stats);
    }

    DEBUG_STATS_PRINT(tasklet_stats, mutex_miscellaneous);

    DPU_TASKLET_STATS_WRITE(&tasklet_stats, (mram_addr_t)(&DPU_TASKLET_STATS_VAR[tasklet_id]));

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
    barrier_id_t barrier = BARRIER_GET(init_barrier);
    perfcounter_t start_time, current_time;
    dpu_compute_time_t accumulate_time;

    if (tasklet_id == 0) {
        accumulate_time = 0ULL;
        mem_reset();
        perfcounter_config(COUNT_CYCLES, true);
        current_time = start_time = perfcounter_get();

        MRAM_INFO_READ(&mram_info, DPU_MRAM_HEAP_POINTER);
        DEBUG_MRAM_INFO_PRINT(&mram_info);

        request_pool_init(&mram_info);
        result_pool_init();

        if (((NB_BYTES_TO_SYMS(mram_info.nbr_len, mram_info.delta) + 2) * 3 * 16) >= 0x10000) {
            printf("cannot run code: symbol length is larger than mulub operation\n");
            halt();
        }

        odpd_init(NB_RUNNING_TASKLETS, NB_BYTES_TO_SYMS(mram_info.nbr_len, mram_info.delta));
    }

    barrier_wait(barrier);

    /* For debugging purpose, one may reduce the number of operating tasklets. */
    if (tasklet_id < NB_RUNNING_TASKLETS) {
        run_align(tasklet_id, &accumulate_time, &current_time);
    }

    barrier_wait(barrier);

    if (tasklet_id == 0) {
        get_time_and_accumulate(&accumulate_time, &current_time);
        accumulate_time = (accumulate_time + current_time) - start_time;
        DPU_COMPUTE_TIME_WRITE(&accumulate_time);
    }

    return 0;
}
