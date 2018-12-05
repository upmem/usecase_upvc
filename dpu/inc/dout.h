/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __DOUT_H__
#define __DOUT_H__

#include "common.h"

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

#define MAX_LOCAL_RESULTS_PER_READ 32
#define LOCAL_RESULTS_PAGE_SIZE (MAX_LOCAL_RESULTS_PER_READ * sizeof(dpu_result_out_t))

/**
 * @brief The results produces by one tasklet when processing one read.

 * @var outs           A local cache of nb_results results.
 * @var nb_results     Total number of results stored.
 * @var mram_base      Base address of the swap area assigned to this tasklet.
 * @var nb_cached_out  Number of results in the local cache.
 * @var nb_page_out    Number of pages of MAX_LOCAL_RESULTS_PER_READ put into the swap area.
 */
typedef struct {
        __attribute__((aligned(8))) dpu_result_out_t outs[MAX_LOCAL_RESULTS_PER_READ];
        unsigned int nb_results;
        mram_addr_t mram_base;
        unsigned int nb_cached_out;
        unsigned int nb_page_out;
} dout_t;

/**
 * @brief Clears the results to start a new request.
 *
 * @param dout  Cleared data output.
 */
void dout_clear(dout_t *dout);

/**
 * @brief Initializes a data out structure for a given tasklet.
 *
 * @param tid the tasklet identifier.
 * @param dout the initialized structure.
 */
void dout_init(unsigned int tid, dout_t *dout);

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
void dout_add(dout_t *dout, uint32_t num, unsigned int score, uint32_t seed_nr, uint32_t seq_nr, dpu_tasklet_stats_t *stats);

/**
 * @brief locates a swap page for a given data out structure.
 *
 * @param dout    Data output.
 * @param pageno  The requested page number.
 */
mram_addr_t dout_swap_page_addr(const dout_t *dout, unsigned int pageno);

#endif /* __DOUT_H__ */
