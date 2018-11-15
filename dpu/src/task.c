/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <defs.h>
#include <stdint.h>
#include <mram.h>
#include <mutex.h>
#include <barrier.h>
#include <alloc.h>
#include <string.h>
#include <mbox.h>

#define NB_BYTES_TO_SYMS(len) ((len) << 2)

#include "debug.h"
#include "odpd.h"
#include "nodp.h"
#include "dout.h"
#include "mutex_def.h"
#include "request_pool.h"
#include "result_pool.h"

#include "common.h"

/**
 * @brief Maximum score allowed.
 */
#define MAX_SCORE (40)

/**
 * @brief The MRAM information, shared between the tasklets.
 *
 * Information are gathered by the first tasklet during boot and used by other tasklets
 * to get the memory topology.
 */
__attribute__((aligned(8))) static mram_info_t mram_info;

/**
 * @brief Fetches a neighbour from the neighbour area.
 *
 * @param base    Offset to the first neighbour of the requested pool within the neighbour area.
 * @param idx     Index of this neighbour in the specified pool.
 * @param buffer  To contain the result, must be the size of a neighbour and coordinates, plus 64 bits alignment.
 * @param stats   To update statistical report.
 *
 * @return A pointer to the beginning of the loaded neighbour.
 */
static uint8_t * load_reference_nbr_and_coords_at(unsigned int base,
                                                  unsigned int idx,
                                                  uint8_t *cache,
                                                  dpu_tasklet_stats_t *stats)
{
        /* The input starts with coordinates (8 bytes), followed by the neighbour. Structure is aligned
         * on 8 bytes boundary.
         */
        unsigned int coords_nbr_len = 8 + ALIGN_DPU(mram_info.nbr_len);
        mram_addr_t coords_nbr_address = (mram_addr_t) (DPU_INPUTS_ADDR + (base + idx) * coords_nbr_len);
        ASSERT_DMA_ADDR(coords_nbr_address, cache, coords_nbr_len);
        ASSERT_DMA_LEN(coords_nbr_len);
        mram_readX(coords_nbr_address, cache, coords_nbr_len);
        stats->mram_data_load += coords_nbr_len;
        stats->mram_load += coords_nbr_len;
        return cache + 8;
}

/**
 * @brief Executes the mapping/align procedure.
 *
 * This function is executed by every tasklet, which picks the requests from the request pool
 * and perform the comparison with the reference genome.
 *
 * @return zero.
 */
static void run_align()
{
        dpu_tasklet_stats_t tasklet_stats = {
                                             .nb_reqs = 0,
                                             .nb_nodp_calls = 0,
                                             .nb_odpd_calls = 0,
                                             .nb_results = 0,
                                             .mram_data_load = 0,
                                             .mram_result_store = 0,
                                             .mram_load = 0,
                                             .mram_store = 0
        };
        mutex_t mutex_miscellaneous = mutex_get(MUTEX_MISCELLANEOUS);
        dout_t dout;
        dpu_request_t request;
        int mini;
        unsigned int nbr_len_aligned = ALIGN_DPU(mram_info.nbr_len);
        DEBUG_STATS_STRUCT;
        DEBUG_RESULTS_VAR;
        DEBUG_REQUESTS_VAR;

        DEBUG_PROCESS_PROFILE(mram_info.nbr_len);
        DEBUG_STATS_CLEAR;

        dout_init(me(), &dout);

        uint8_t *current_read_nbr = mem_alloc(nbr_len_aligned);
        memset(current_read_nbr, 0, nbr_len_aligned);

        /* Create a local cache to get the reference reads, with few more bytes for address
         * alignment.
         */
        uint8_t *cached_coords_and_nbr = mem_alloc_dma(nbr_len_aligned + 2 * sizeof(uint32_t) + 16);
        memset(cached_coords_and_nbr, 0, nbr_len_aligned + 2 * sizeof(uint32_t) + 16);

        while (request_pool_next(&request, current_read_nbr, &tasklet_stats, &mram_info)) {

                DEBUG_REQUESTS_PRINT(request, mram_info.nbr_len);
                DEBUG_RESULTS_INCR;
                DEBUG_STATS_INC_NB_READS;

                tasklet_stats.nb_reqs++;

                mini = MAX_SCORE;
                dout_clear(&dout);

                for (unsigned int idx = 0; idx < request.count; idx++) {
                        int score, score_nodp, score_odpd = -1;
                        uint8_t *ref_nbr = load_reference_nbr_and_coords_at(request.offset, idx, cached_coords_and_nbr,
                                                                            &tasklet_stats);

                        DEBUG_REQUESTS_PRINT_REF(ref_nbr, mram_info.nbr_len);
                        DEBUG_STATS_INC_NB_REFS;

                        score = score_nodp = noDP(current_read_nbr, ref_nbr, mram_info.nbr_len, mini);
                        tasklet_stats.nb_nodp_calls++;

                        if (score_nodp == -1) {
                                score_odpd = score = odpd(current_read_nbr,
                                                              ref_nbr,
                                                              mini,
                                                              NB_BYTES_TO_SYMS(mram_info.nbr_len),
                                                              me());

                                tasklet_stats.nb_odpd_calls++;
                                DEBUG_STATS_INC_NB_ODPD_CALLS;
                        }
                        DEBUG_PROCESS_SCORES(mram_info.nbr_len,
                                             current_read_nbr,
                                             ref_nbr,
                                             mini,
                                             score_nodp,
                                             score_odpd,
                                             mutex_miscellaneous);

                        if (score <= mini) {
                                if (score < mini) {
                                        mini = score;
                                        /* Get rid of previous results, we found a better one */
                                        dout_clear(&dout);
                                }
                                if (dout.nb_results < MAX_RESULTS_PER_READ) {
                                        dout_add(&dout, request.num, (unsigned int) score,
                                                 ((uint32_t *) cached_coords_and_nbr)[0],
                                                 ((uint32_t *) cached_coords_and_nbr)[1],
                                                 &tasklet_stats
                                                 );
                                        DEBUG_RESULTS_PRINT(me(),
                                                            idx,
                                                            request.num,
                                                            request.offset,
                                                            ((uint32_t *) cached_coords_and_nbr)[0],
                                                            ((uint32_t *) cached_coords_and_nbr)[1],
                                                            score);
                                } else {
                                        printf("WARNING! too many results for request!\n");
                                        /* Trigger a fault, since this should never happen. */
                                        halt();
                                }
                        }
                }
                if (dout.nb_results != 0) {
                        tasklet_stats.nb_results += dout.nb_results;
                        result_pool_write(&dout, &tasklet_stats);
                }
        }

        DEBUG_STATS_PRINT;
        DEBUG_STATS_MRAM_PRINT(request_pool_get_stats_load(), result_pool_get_stats_write(), mutex_miscellaneous);

        mbox_send(&tasklet_stats, sizeof(tasklet_stats));
}

/**
 * @brief Common main: gather the seed mapping for this DPU, then start processing the requests on every tasklet.
 *
 * @return 0 if everything worked fine
 */
int main()
{
        barrier_t barrier = barrier_get(0);
        if (me() == 0) {
                MRAM_INFO_READ(MRAM_INFO_ADDR, &mram_info);
                DEBUG_MRAM_INFO_PRINT(&mram_info);

                request_pool_init(&mram_info);
                result_pool_init();

                if (((NB_BYTES_TO_SYMS(mram_info.nbr_len) + 1) * 24) >= 4096) {
                        printf("cannot run code: symbol length is larger than mulub operation\n");
                        halt();
                }

                odpd_init(NB_RUNNING_TASKLETS, NB_BYTES_TO_SYMS(mram_info.nbr_len));
        }
        barrier_wait(barrier);
        /* For debugging purpose, one may reduce the number of operating tasklets. */
        if (me() < NB_RUNNING_TASKLETS) {
                run_align();
        }

        return 0;
}
