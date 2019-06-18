/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <alloc.h>
#include <defs.h>
#include <mram.h>

#include "dout.h"
#include "debug.h"
#include "stats.h"

__mram dpu_result_out_t m_dpu_swap_result[NB_TASKLET_PER_DPU * MAX_RESULTS_PER_READ];

void dout_clear(dout_t *dout)
{
    dout->nb_results = 0;
    dout->nb_page_out = 0;
    dout->nb_cached_out = 0;
}

void dout_init(unsigned int tid, dout_t *dout)
{
    dout->mram_base = (mram_addr_t)&m_dpu_swap_result[tid * MAX_RESULTS_PER_READ];
    dout_clear(dout);
}

void dout_add(dout_t *dout, uint32_t num, unsigned int score, uint32_t seed_nr, uint32_t seq_nr, dpu_tasklet_stats_t *stats)
{
    dpu_result_out_t *new_out;
    if (dout->nb_cached_out == MAX_LOCAL_RESULTS_PER_READ) {
        mram_addr_t swap_addr;

        /* Local cache is full, copy into the swap area.
         * This shall never happen, but let's verify that we remain inside our assigned swap area. */
        if (dout->nb_page_out >= MAX_RESULTS_PER_READ / MAX_LOCAL_RESULTS_PER_READ) {
            printf("WARNING! too many swapped pages!\n");
            halt();
        }
        swap_addr = dout_swap_page_addr(dout, dout->nb_page_out);
        LOCAL_RESULTS_PAGE_WRITE(dout->outs, swap_addr);
        STATS_INCR_STORE(stats, LOCAL_RESULTS_PAGE_SIZE);
        dout->nb_cached_out = 0;
        dout->nb_page_out++;
    }

    new_out = &dout->outs[dout->nb_cached_out];
    new_out->num = num;
    new_out->score = score;
    new_out->coord.seed_nr = seed_nr;
    new_out->coord.seq_nr = seq_nr;

    dout->nb_cached_out++;
    dout->nb_results++;
}

mram_addr_t dout_swap_page_addr(const dout_t *dout, unsigned int pageno)
{
    return (mram_addr_t)(dout->mram_base + (pageno * LOCAL_RESULTS_PAGE_SIZE));
}
