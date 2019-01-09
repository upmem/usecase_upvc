/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <mram.h>
#include <alloc.h>
#include <defs.h>
#include <mutex.h>
#include <string.h>

#include "debug.h"
#include "mutex_def.h"
#include "dout.h"
#include "stats.h"

#include "common.h"

/**
 * @brief Common structure to write back results.
 *
 * This output FIFO is shared by the tasklets to write back results to host, thus is protected by a critical section.
 *
 * @var cache      Local cache to perform memory transfers.
 * @var mutex      Critical section that protects the pool.
 * @var wridx      Index of the current output in the FIFO.
 * @var cur_write  Where to write in MRAM.
 */
typedef struct {
        uint8_t cache[LOCAL_RESULTS_PAGE_SIZE];
        mutex_t mutex;
        unsigned int wridx;
        mram_addr_t cur_write;
} result_pool_t;

/**
 * @brief The result pool shared by tasklets.
 */
__attribute__((aligned(8))) static result_pool_t result_pool;

void result_pool_init()
{
        result_pool.mutex = mutex_get(MUTEX_RESULT_POOL);
        result_pool.wridx = 0;
        result_pool.cur_write = DPU_RESULT_ADDR;
}

void result_pool_write(const dout_t *results, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
        unsigned int each_result = 0;
        unsigned int pageno;

        mutex_lock(result_pool.mutex);

        /* Read back and write the swapped results */
        for (pageno = 0; pageno < results->nb_page_out; pageno++) {
                mram_addr_t source_addr = dout_swap_page_addr(results, pageno);

                if (result_pool.wridx + MAX_LOCAL_RESULTS_PER_READ >= (MAX_DPU_RESULTS - 1)) {
                        printf("WARNING! too many result in DPU! (from swap)\n");
                        halt();
                }

                ASSERT_DMA_ADDR(source_addr, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
                STATS_INCR_LOAD(stats, LOCAL_RESULTS_PAGE_SIZE);
                LOCAL_RESULTS_PAGE_READ(source_addr, result_pool.cache);

                ASSERT_DMA_ADDR(result_pool.cur_write, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
                STATS_INCR_STORE(stats, LOCAL_RESULTS_PAGE_SIZE);
                STATS_INCR_STORE_RESULT(stats, LOCAL_RESULTS_PAGE_SIZE);
                LOCAL_RESULTS_PAGE_WRITE(result_pool.cache, result_pool.cur_write);

                result_pool.wridx += MAX_LOCAL_RESULTS_PER_READ;
                result_pool.cur_write += LOCAL_RESULTS_PAGE_SIZE;
        }

        while ((each_result < results->nb_cached_out) && (result_pool.wridx < (MAX_DPU_RESULTS - 1))) {
                /* Ensure that the size of a result out structure is two longs. */
                ASSERT_DMA_ADDR(result_pool.cur_write, &(results->outs[each_result]), sizeof(dpu_result_out_t));
                STATS_INCR_STORE(stats, sizeof(dpu_result_out_t));
                STATS_INCR_STORE_RESULT(stats, sizeof(dpu_result_out_t));
                DPU_RESULT_WRITE((void *)&(results->outs[each_result]), result_pool.cur_write);

                result_pool.wridx++;
                result_pool.cur_write += sizeof(dpu_result_out_t);
                each_result++;
        }
        if (result_pool.wridx >= (MAX_DPU_RESULTS - 1)) {
                printf("WARNING! too many result in DPU! (from local)\n");
                halt();
        }

        mutex_unlock(result_pool.mutex);
}

void result_pool_finish(STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
        __attribute__((aligned(8))) static const dpu_result_out_t end_of_results =
                {
                 .num = (unsigned int) -1,
                 .score = (unsigned int) -1,
                 .coord.seq_nr = 0,
                 .coord.seed_nr = 0
                };
        /* Note: will fill in the result pool until MAX_DPU_RESULTS -1, to be sure that the very last result
         * has a num equal to -1.
         */
        mutex_lock(result_pool.mutex);
        /* Mark the end of result data, do not increment the indexes, so that the next one restarts from this
         * point.
         */
        DPU_RESULT_WRITE((void*)&end_of_results, result_pool.cur_write);
        STATS_INCR_STORE(stats, sizeof(dpu_result_out_t));
        STATS_INCR_STORE_RESULT(stats, sizeof(dpu_result_out_t));

        mutex_unlock(result_pool.mutex);
}
