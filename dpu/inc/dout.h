/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DOUT_H__
#define __INTEGRATION_DOUT_H__

#include <mram.h>
#include "debug.h"

/**
 * @brief Management of the output data produced by tasklets.
 *
 * One tasklet can produce up to MAX_RESULTS_PER_READ for one request. Empirically, we observe that the
 * value can be up to more than 500 results. Given that a result is 16 bytes (type result_out_t),
 * this represents an amount of almost 512x16x16=128KB, which does not give enough space
 * to the rest of the application.
 *
 * As a consequence, the production of results is relayed by a swapping system, to store data into MRAM.
 * The MRAM holds a buffer before the output area, which can contain up to 1024 results, representing
 * 1024x16x16=256KB. This "swap area" sits in the 256KB preceding the output area,
 * i.e. 64MB-16xMAX_DPU_RESULTS-256K.
 *
 * The "data out" (dout) module manages both the local caching of results and the swapping.
 */

/* See mdpu.h */
#define MAX_DPU_RESULTS (1 << 16)
#define MRAM_SIZE (64 << 20)

/**
 * @brief The maximum number of results produced locally before swapping to MRAM.
 */
#define MAX_LOCAL_RESULTS_PER_READ 32

/**
 * @brief One "page" of results: 1 page is cached, the olders are swapped in MRAM.
 */
#define LOCAL_RESULTS_PAGE_SIZE (MAX_LOCAL_RESULTS_PER_READ * sizeof(result_out_t))

/**
 * @brief The maximum number of results produced for one read.
 */
#define MAX_RESULTS_PER_READ 1024

/**
 * @brief One result produced for one read
 *
 * Note: please be careful to keep the size of this structure as a multiple of 64 bits, to optimize MRAM
 * transfers.
 *
 * See dpus.h/dpu_result_out_t for consistency of definitions.
 */
typedef struct {
        unsigned int num;
        uint32_t seq_nr;
        uint32_t seed_nr;
        unsigned int score;
} result_out_t;

/**
 * @brief stats reported by every tasklet
 *
 * Note: the structure size must be consistent with tasklet mailbox sizes. Please review build_configuration.script
 * if changing this.
 *
 * See dpus.h/dpu_stats_out_t for consistency of definitions.
 */
typedef struct {
        uint32_t nb_reqs;
        uint32_t nb_nodp_calls;
        uint32_t nb_odpd_calls;
        uint32_t nb_results;
        uint32_t mram_data_load;
        uint32_t mram_result_store;
        uint32_t mram_load;
        uint32_t mram_store;
} dpu_tasklet_stats_t;

/**
 * @brief The results produces by one tasklet when processing one read.

 * @var nb_results     Total number of results stored.
 * @var outs           A local cache of nb_results results.
 * @var mram_base      Base address of the swap area assigned to this tasklet.
 * @var nb_cached_out  Number of results in the local cache.
 * @var nb_page_out    Number of pages of MAX_LOCAL_RESULTS_PER_READ put into the swap area.
 */
typedef struct {
        unsigned int nb_results;
        result_out_t *outs;
        mram_addr_t mram_base;
        unsigned int nb_cached_out;
        unsigned int nb_page_out;
} dout_t;

/**
 * @brief Clears the results to start a new request.
 *
 * @param dout  Cleared data output.
 */
static void dout_clear(dout_t *dout)
{
        dout->nb_results = 0;
        dout->nb_page_out = 0;
        dout->nb_cached_out = 0;
}

/**
 * @brief Initializes a data out structure for a given tasklet.
 *
 * @param tid the tasklet identifier.
 * @param dout the initialized structure.
 */
static void dout_init(unsigned int tid, dout_t *dout)
{
        dout->mram_base = MRAM_SIZE -
                (MAX_DPU_RESULTS * sizeof(result_out_t)) -
                (16 * MAX_RESULTS_PER_READ * sizeof(result_out_t)) +
                (tid * MAX_RESULTS_PER_READ * sizeof(result_out_t));
        dout->outs = mem_alloc_dma(LOCAL_RESULTS_PAGE_SIZE);
        dout_clear(dout);
}

/**
 * @brief Tecords a new result
 *
 * @param dout     Data output.
 * @param num      Request number.
 * @param score    Recorded scored.
 * @param seed_nr  Recorded seed number.
 * @param seq_nr   Recorded sequence number.
 * @param stats    To update statistical report.
 */
static void dout_add(dout_t *dout, uint32_t num, unsigned int score, uint32_t seed_nr, uint32_t seq_nr, dpu_tasklet_stats_t *stats)
{
        result_out_t *new_out;
        if (dout->nb_cached_out == MAX_LOCAL_RESULTS_PER_READ) {
                mram_addr_t swap_addr;

                /* Local cache is full, copy into the swap area.
                 * This shall never happen, but let's verify that we remain inside our assigned swap area. */
                if (dout->nb_page_out >= MAX_RESULTS_PER_READ / MAX_LOCAL_RESULTS_PER_READ) {
                        printf("WARNING! too many swapped pages!\n");
                        halt();
                }
                swap_addr = dout->mram_base + dout->nb_page_out * LOCAL_RESULTS_PAGE_SIZE;
                mram_writeX(dout->outs, swap_addr, LOCAL_RESULTS_PAGE_SIZE);
                stats->mram_store += LOCAL_RESULTS_PAGE_SIZE;
                dout->nb_cached_out = 0;
                dout->nb_page_out++;
        }

        new_out = dout->outs + dout->nb_cached_out;
        new_out->num = num;
        new_out->score = score;
        new_out->seed_nr = seed_nr;
        new_out->seq_nr = seq_nr;

        dout->nb_cached_out++;
        dout->nb_results++;
}

/**
 * @brief locates a swap page for a given data out structure.
 *
 * @param dout    Data output.
 * @param pageno  The requested page number.
 */
static mram_addr_t dout_swap_page_addr(const dout_t *dout, unsigned int pageno)
{
        return (mram_addr_t) (dout->mram_base + (pageno * LOCAL_RESULTS_PAGE_SIZE));
}

#endif /* __INTEGRATION_DOUT_H__ */
