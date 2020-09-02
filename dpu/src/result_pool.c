/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <alloc.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <string.h>

#include "debug.h"
#include "dout.h"
#include "stats.h"

#include "common.h"

/**
 * @brief Common structure to write back results.
 *
 * This output FIFO is shared by the tasklets to write back results to host, thus is protected by a critical section.
 *
 * @var cache      Local cache to perform memory transfers.
 * @var cur_write  Where to write in MRAM.
 */
typedef struct {
    uint8_t cache[LOCAL_RESULTS_PAGE_SIZE];
    uintptr_t cur_write;
} result_pool_t;

__host nb_result_t DPU_NB_RESULT_VAR0;
__host nb_result_t DPU_NB_RESULT_VAR1;
__host nb_result_t DPU_NB_RESULT_VAR2;
__host nb_result_t DPU_NB_RESULT_VAR3;
nb_result_t *DPU_NB_RESULT_VAR;

/**
 * @brief The result pool shared by tasklets.
 */
__dma_aligned static result_pool_t result_pool;
MUTEX_INIT(result_pool_mutex);

/**
 * @brief The buffer of result in mram.
 */
__mram_noinit dpu_result_out_t DPU_RESULT_VAR0[MAX_DPU_RESULTS];
__mram_noinit dpu_result_out_t DPU_RESULT_VAR1[MAX_DPU_RESULTS];
__mram_noinit dpu_result_out_t DPU_RESULT_VAR2[MAX_DPU_RESULTS];
__mram_noinit dpu_result_out_t DPU_RESULT_VAR3[MAX_DPU_RESULTS];

__mram_ptr dpu_result_out_t *DPU_RESULT_VAR;
__host uint32_t dpu_result_var_idx;

void result_pool_init()
{
    if (dpu_result_var_idx%4 ==0) {
        DPU_NB_RESULT_VAR = &DPU_NB_RESULT_VAR0;
        DPU_RESULT_VAR = DPU_RESULT_VAR0;
    } else if (dpu_result_var_idx%4 == 1){
        DPU_NB_RESULT_VAR = &DPU_NB_RESULT_VAR1;
        DPU_RESULT_VAR = DPU_RESULT_VAR1;
    } else if (dpu_result_var_idx%4 == 2){
        DPU_NB_RESULT_VAR = &DPU_NB_RESULT_VAR2;
        DPU_RESULT_VAR = DPU_RESULT_VAR2;
    } else if (dpu_result_var_idx%4 == 3){
        DPU_NB_RESULT_VAR = &DPU_NB_RESULT_VAR3;
        DPU_RESULT_VAR = DPU_RESULT_VAR3;
    }
    *DPU_NB_RESULT_VAR = 0;
    result_pool.cur_write = (uintptr_t)DPU_RESULT_VAR;
}

void result_pool_write(const dout_t *results, STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    unsigned int each_result = 0;
    unsigned int pageno;

    mutex_lock(result_pool_mutex);

    /* Read back and write the swapped results */
    for (pageno = 0; pageno < results->nb_page_out; pageno++) {
        __mram_ptr void *source_addr = dout_swap_page_addr(results, pageno);

        if (*DPU_NB_RESULT_VAR + MAX_LOCAL_RESULTS_PER_READ >= (MAX_DPU_RESULTS - 1)) {
            printf("WARNING! too many result in DPU! (from swap)\n");
            halt();
        }

        ASSERT_DMA_ADDR(source_addr, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
        STATS_INCR_LOAD(stats, LOCAL_RESULTS_PAGE_SIZE);
        mram_read(source_addr, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);

        ASSERT_DMA_ADDR(result_pool.cur_write, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
        STATS_INCR_STORE(stats, LOCAL_RESULTS_PAGE_SIZE);
        STATS_INCR_STORE_RESULT(stats, LOCAL_RESULTS_PAGE_SIZE);
        mram_write(result_pool.cache, (__mram_ptr void *)result_pool.cur_write, LOCAL_RESULTS_PAGE_SIZE);

        *DPU_NB_RESULT_VAR += MAX_LOCAL_RESULTS_PER_READ;
        result_pool.cur_write += LOCAL_RESULTS_PAGE_SIZE;
    }

    while ((each_result < results->nb_cached_out) && (*DPU_NB_RESULT_VAR < (MAX_DPU_RESULTS - 1))) {
        /* Ensure that the size of a result out structure is two longs. */
        ASSERT_DMA_ADDR(result_pool.cur_write, &(results->outs[each_result]), sizeof(dpu_result_out_t));
        STATS_INCR_STORE(stats, sizeof(dpu_result_out_t));
        STATS_INCR_STORE_RESULT(stats, sizeof(dpu_result_out_t));
        mram_write((void *)&(results->outs[each_result]), (__mram_ptr void *)result_pool.cur_write, sizeof(dpu_result_out_t));

        *DPU_NB_RESULT_VAR += 1;
        result_pool.cur_write += sizeof(dpu_result_out_t);
        each_result++;
    }
    if (*DPU_NB_RESULT_VAR >= (MAX_DPU_RESULTS - 1)) {
        printf("WARNING! too many result in DPU! (from local)\n");
        halt();
    }

    mutex_unlock(result_pool_mutex);
}

void result_pool_finish(STATS_ATTRIBUTE dpu_tasklet_stats_t *stats)
{
    __dma_aligned static const dpu_result_out_t end_of_results
        = { .num = (unsigned int)-1, .score = (unsigned int)-1, .coord.seq_nr = 0, .coord.seed_nr = 0 };
    /* Note: will fill in the result pool until MAX_DPU_RESULTS -1, to be sure that the very last result
     * has a num equal to -1.
     */
    mutex_lock(result_pool_mutex);
    /* Mark the end of result data, do not increment the indexes, so that the next one restarts from this
     * point.
     */
    mram_write((void *)&end_of_results, (__mram_ptr void *)result_pool.cur_write, sizeof(dpu_result_out_t));
    STATS_INCR_STORE(stats, sizeof(dpu_result_out_t));
    STATS_INCR_STORE_RESULT(stats, sizeof(dpu_result_out_t));

    mutex_unlock(result_pool_mutex);
}

void check_end_mark() {
    __dma_aligned dpu_result_out_t res;
    mram_read((__mram_ptr void *)((64 << 20) - sizeof(dpu_result_out_t)), (void *)&res, sizeof(dpu_result_out_t));
    mram_read((__mram_ptr void *)&DPU_RESULT_VAR0[DPU_NB_RESULT_VAR0], (void *)&res, sizeof(dpu_result_out_t));
    if (res.num != -1 || res.score != -1 || res.coord.seq_nr != 0 || res.coord.seed_nr != 0) {
        halt();
    }
    mram_read((__mram_ptr void *)((64 << 20) - sizeof(dpu_result_out_t)), (void *)&res, sizeof(dpu_result_out_t));
    mram_read((__mram_ptr void *)&DPU_RESULT_VAR1[DPU_NB_RESULT_VAR1], (void *)&res, sizeof(dpu_result_out_t));
    if (res.num != -1 || res.score != -1 || res.coord.seq_nr != 0 || res.coord.seed_nr != 0) {
        halt();
    }
    mram_read((__mram_ptr void *)((64 << 20) - sizeof(dpu_result_out_t)), (void *)&res, sizeof(dpu_result_out_t));
    mram_read((__mram_ptr void *)&DPU_RESULT_VAR2[DPU_NB_RESULT_VAR2], (void *)&res, sizeof(dpu_result_out_t));
    if (res.num != -1 || res.score != -1 || res.coord.seq_nr != 0 || res.coord.seed_nr != 0) {
        halt();
    }
    mram_read((__mram_ptr void *)((64 << 20) - sizeof(dpu_result_out_t)), (void *)&res, sizeof(dpu_result_out_t));
    mram_read((__mram_ptr void *)&DPU_RESULT_VAR3[DPU_NB_RESULT_VAR3], (void *)&res, sizeof(dpu_result_out_t));
    if (res.num != -1 || res.score != -1 || res.coord.seq_nr != 0 || res.coord.seed_nr != 0) {
        halt();
    }
}
