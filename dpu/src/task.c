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
#include <request.h>
#include <mbox.h>
#include <dout.h>

#include "debug.h"
#include "odpd.h"
#include "nodp.h"
#include "dout.h"

#define MUTEX_REQUEST_POOL  0
#define MUTEX_RESULT_POOL   1
#define MUTEX_MISCELLANEOUS 2

/**
 * @brief Maximum score allowed.
 */
#define MAX_SCORE (40)

/**
 * See mdpu.h. Notice that the host ensures that every area is aligned on 64 bits.
 */
typedef struct {
        uint32_t io_offs;
        uint32_t usage;
        uint32_t nb_nbr;
        uint32_t nbr_offs;
        uint32_t nbr_len;
        uint32_t unused;
} mram_info_t;

/**
 * @brief The MRAM information, shared between the tasklets.
 *
 * Information are gathered by the first tasklet during boot and used by other tasklets
 * to get the memory topology.
 */
static mram_info_t *mram_info;

#ifdef DEBUG_MRAM_INFO

static void print_mram_info() {
        printf("MRAM usage       = %u\n", mram_info->usage);
        printf("MRAM nr neigbors = %u\n", mram_info->nb_nbr);
        printf("MRAM nbr len     = %u\n", mram_info->nbr_len);
        printf("MRAM nbr@0x%x\n", mram_info->nbr_offs);
        printf("MRAM I/O@0x%x\n", mram_info->io_offs);
        printf("MRAM magic=0x%x\n", mram_info->unused);
}

#else
#define print_mram_info() do {} while(0)
#endif

/**
 * @brief Loads the base information used by this application to work with MRAM data.
 *
 * The shared variable mram_info is updated with the header information.
 */
static void load_mram_info() {
        mram_info = mem_alloc_dma(sizeof(mram_info_t));
        mram_read24((mram_addr_t) 0, mram_info);
}

/**
 * @brief Common structure to consume requests.
 *
 * Requests belong to a FIFO, from which each tasklet picks the reads. The FIFO is protected by a critical
 * section.
 *
 * @var mutex       Critical section that protects the pool.
 * @var nb_reads    The number of reads in the request pool.
 * @var rdidx       Index of the first unread read in the request pool.
 * @var cur_read    Address of the first read to be processed in MRAM.
 * @var cache       An internal cache to fetch requests from MRAM.
 * @var cache_size  The cache size
 */
typedef struct {
        mutex_t mutex;
        unsigned int nb_reads;
        unsigned int rdidx;
        mram_addr_t cur_read;
        unsigned int cache_size;
        uint8_t *cache;
        unsigned int stats_load;
} request_pool_t;

/**
 * @brief Common request pool, shared by every tasklet.
 */
request_pool_t request_pool;

/**
 * @brief Common structure to write back results.
 *
 * This output FIFO is shared by the tasklets to write back results to host, thus is protected by a critical section.
 *
 * @var mutex      Critical section that protects the pool.
 * @var wridx      Index of the current output in the FIFO.
 * @var cur_write  Where to write in MRAM.
 * @var cache      Local cache to perform memory transfers.
 */
typedef struct {
        mutex_t mutex;
        unsigned int wridx;
        mram_addr_t cur_write;
        uint8_t *cache;
        unsigned int stats_write;
} result_pool_t;

/**
 * @brief The result pool shared by tasklets.
 */
result_pool_t result_pool;

/**
 * @brief Prints the current state of the request pool.
 */
#ifdef DEBUG_REQUESTS

static void request_pool_print() {
        printf("R> nb_reads = %u\n", request_pool.nb_reads);
        printf("R> rdidx    = %u\n", request_pool.rdidx);
        printf("R> cur_read = %x\n", (uint32_t) (request_pool.cur_read));
}

#else
#define request_pool_print() do {} while(0)
#endif // DEBUG_REQUESTS

/**
 * @return The transfer size needed to load a neighbour as a byte stream.
 */
static inline unsigned int nbr_len_as_longs() {
        return (mram_info->nbr_len + 7) & ~7;
}

/**
 * @brief Initializes the request pool.
 * @param io_offs  The beginning of the I/O area in MRAM.
 */
static void request_pool_init(mram_addr_t io_offs) {
        /* The I/O area header, as written by the host in dpu_try_write_dispatch_into_mram. */
        struct {
                uint32_t nb_reads;
                uint32_t magic;
        } *io_data = mem_alloc_dma(2 * sizeof(uint32_t));

        request_pool.mutex = mutex_get(MUTEX_REQUEST_POOL);

        mram_read8((mram_addr_t) (io_offs), io_data);
        request_pool.nb_reads = io_data->nb_reads;
        request_pool.rdidx = 0;
        request_pool.cur_read = (mram_addr_t) (io_offs + sizeof(*io_data));
        /* Every request is aligned on 64 bits. Consequently, we must reserve a cache
         * able to contain a neighbour and a request header, aligned on a multiple of
         * 8 bytes.
         */
        request_pool.cache_size = (SIZEOF_REQUEST_HEADER + mram_info->nbr_len + 7) & ~7;
        request_pool.cache = (uint8_t *) mem_alloc_dma(request_pool.cache_size);
        request_pool.stats_load = 0;
        request_pool_print();
}

/**
 * @brief Gets the next read from the request pool, if any.
 *
 * @param request  Output Request filled with the next request, if found.
 * @param nbr      Ouput neighbour filled with the requested neighbour, if found.
 * @param stats    To update statistical reports.
 *
 * @return True if a new request was fetched, false if the FIFO is empty.
 */
static bool request_pool_next(request_t *request, uint8_t *nbr, dpu_tasklet_stats_t *stats) {
        mutex_lock(request_pool.mutex);
        if (request_pool.rdidx == request_pool.nb_reads) {
                mutex_unlock(request_pool.mutex);
                return false;
        }

        /* Fetch next request into cache */
        assert_dma_addr(request_pool.cur_read, request_pool.cache, request_pool.cache_size);
        assert_dma_len(request_pool.cache_size);
        mram_readX(request_pool.cur_read, (void *) request_pool.cache, request_pool.cache_size);
        stats->mram_load += request_pool.cache_size;
        stats->mram_data_load += request_pool.cache_size;
        request_pool.stats_load += request_pool.cache_size;

        memcpy(request, request_pool.cache, sizeof(request_t));
        memcpy(nbr, request_pool.cache + SIZEOF_REQUEST_HEADER, mram_info->nbr_len);

        /* Point to next request */
        request_pool.rdidx++;
        request_pool.cur_read += request_pool.cache_size;
        mutex_unlock(request_pool.mutex);

        return true;
}

/**
 * @brief Initializes the result pool.
 *
 * @param out_offset  Start address of this pool in MRAM.
 */
static void result_pool_init(mram_addr_t out_offset) {
        /* Preliminary security: if the size of result_out changed in the code, result_pool_write does not work
         * anymore.
         */
        if (sizeof(result_out_t) != 16) {
                printf("!!!! the size of result_out changed, please update result_pool functions\n");
                halt();
        }

        result_pool.mutex = mutex_get(MUTEX_RESULT_POOL);
        result_pool.wridx = 0;
        result_pool.cur_write = out_offset;
        /* Will write the results by pages */
        result_pool.cache = (uint8_t *) mem_alloc_dma(LOCAL_RESULTS_PAGE_SIZE);
        result_pool.stats_write = 0;
}

/**
 * @brief Pushes a bunch of results.
 *
 * @param results  The list of outputs.
 * @param stats    To update statistical report.
 */
static void result_pool_write(const dout_t *results, dpu_tasklet_stats_t *stats) {
        mutex_lock(result_pool.mutex);
        unsigned int i;
        unsigned int pageno;
        static result_out_t end_of_results = {
                                              .num = (unsigned int) -1,
                                              .score = (unsigned int) -1,
                                              .seq_nr = 0,
                                              .seed_nr = 0
        };
        /* Note: will fill in the result pool until MAX_DPU_RESULTS -1, to be sure that the very last result
         * has a num equal to -1.
         */

        /* Read back and write the swapped results */
        for (pageno = 0; pageno < results->nb_page_out; pageno++) {
                mram_addr_t source_addr = dout_swap_page_addr(results, pageno);
                assert_dma_addr(source_addr, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
                mram_readX(source_addr, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
                stats->mram_load += LOCAL_RESULTS_PAGE_SIZE;
                if (result_pool.wridx + MAX_LOCAL_RESULTS_PER_READ < (MAX_DPU_RESULTS - 1)) {
                        assert_dma_addr(result_pool.cur_write, result_pool.cache, LOCAL_RESULTS_PAGE_SIZE);
                        mram_writeX(result_pool.cache, result_pool.cur_write, LOCAL_RESULTS_PAGE_SIZE);
                        stats->mram_store += LOCAL_RESULTS_PAGE_SIZE;
                        stats->mram_result_store += LOCAL_RESULTS_PAGE_SIZE;
                        result_pool.wridx += MAX_LOCAL_RESULTS_PER_READ;
                        result_pool.cur_write += LOCAL_RESULTS_PAGE_SIZE;
                        result_pool.stats_write += LOCAL_RESULTS_PAGE_SIZE;
                }
        }

        for (i = 0; (result_pool.wridx < (MAX_DPU_RESULTS - 1)) && (i < results->nb_cached_out); i++) {
                (void) memcpy(result_pool.cache, &(results->outs[i]), sizeof(result_out_t));
                /* Ensure that the size of a result out structure is two longs. */
                assert_dma_addr(result_pool.cur_write, result_pool.cache, 16);
                mram_write16(result_pool.cache, result_pool.cur_write);
                stats->mram_store += 16;
                stats->mram_result_store += 16;
                result_pool.wridx++;
                result_pool.cur_write += sizeof(result_out_t);
                result_pool.stats_write += 16;
        }

        /* Mark the end of result data, do not increment the indexes, so that the next one restarts from this
         * point.
         */
        mram_write16(&end_of_results, result_pool.cur_write);
        result_pool.stats_write += 16;
        mutex_unlock(result_pool.mutex);
}

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
static uint8_t *
load_reference_nbr_and_coords_at(unsigned int base, unsigned int idx, uint8_t *cache, dpu_tasklet_stats_t *stats) {
        /* The input starts with coordinates (8 bytes), followed by the neighbour. Structure is aligned
         * on 8 bytes boundary.
         */
        unsigned int coords_nbr_len = 8 + nbr_len_as_longs();
        mram_addr_t coords_nbr_address = (mram_addr_t) (mram_info->nbr_offs + (base + idx) * coords_nbr_len);
        assert_dma_addr(coords_nbr_address, cache, coords_nbr_len);
        assert_dma_len(coords_nbr_len);
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
static int run_align() {
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

#ifdef DEBUG_STATS
        struct stats stats;
#endif
        dout_t dout;
        dout_init(me(), &dout);

        request_t request;
        int mini, score;

        stat_clear(&stats);
#ifdef DEBUG_RESULTS
        /* This debug counter grows the same way the reference code does, so that we can
         * easily compare results via this unique reference number.
         */
        int debug_nr = -1;
#endif /* DEBUG_RESULTS */

        uint8_t *current_read_nbr = mem_alloc(nbr_len_as_longs());
        memset(current_read_nbr, 0, nbr_len_as_longs());

#ifdef DEBUG_REQUESTS
        /* To print a sequence in plain text: unwrap bytes to symbols */
        uint8_t *syms = mem_alloc(nb_bytes_to_syms(nbr_len_as_longs()));
#endif /* DEBUG_REQUESTS */

        /* Create a local cache to get the reference reads, with few more bytes for address
         * alignment.
         */
        uint8_t *cached_coords_and_nbr = mem_alloc_dma(nbr_len_as_longs() + 2 * sizeof(uint32_t) + 16);
        memset(cached_coords_and_nbr, 0, nbr_len_as_longs() + 2 * sizeof(uint32_t) + 16);
        trace_pattern_profile(mram_info->nbr_len);

        while (request_pool_next(&request, current_read_nbr, &tasklet_stats)) {
#ifdef DEBUG_REQUESTS
                decode_nbr((uint8_t *) &(request.nbr), syms, mram_info->nbr_len);
                debug_request(&request, mram_info->nbr_len, syms);
#endif /* DEBUG_REQUESTS */
#ifdef DEBUG_RESULTS
                debug_nr++;
#endif /* DEBUG_RESULTS */

                stat_inc_nb_reads(&stats);
                tasklet_stats.nb_reqs++;

                mini = MAX_SCORE;
                dout_clear(&dout);

                unsigned int idx;
                for (idx = 0; idx < request.count; idx++) {
                        uint8_t *ref_nbr = load_reference_nbr_and_coords_at(request.offset, idx, cached_coords_and_nbr,
                                                                            &tasklet_stats);

#ifdef DEBUG_REQUESTS
                        decode_nbr(ref_nbr, syms, mram_info->nbr_len);
                        debug_reference(mram_info->nbr_len, syms);
#endif /* DEBUG_REQUESTS */


                        int score_nodp = noDP(current_read_nbr, ref_nbr, mram_info->nbr_len, mini);
                        tasklet_stats.nb_nodp_calls++;
                        stat_inc_nb_refs(&stats);

#ifdef DEBUG_PROCESS
                        uint32_t seed_nr = ((uint32_t *) cached_coords_and_nbr)[0];
                        uint32_t seq_nr = ((uint32_t *) cached_coords_and_nbr)[1];
                        printf("%u.%u noDP score = %u\n", seq_nr, seed_nr, score);
#endif

                        if (score_nodp == -1) {
                                stat_inc_nb_odpd_calls(&stats);
#ifndef OPT_S
                                int score_odpd = odpd(current_read_nbr, ref_nbr, mini, nb_bytes_to_syms(mram_info->nbr_len), me());
#else
                                int score_odpd = odpdx(current_read_nbr, ref_nbr, mini, nb_bytes_to_syms(mram_info->nbr_len), me());
#endif /* OPT_S */
                                tasklet_stats.nb_odpd_calls++;
                                score = score_odpd;
                                trace_pattern_and_scores(mram_info->nbr_len, current_read_nbr, ref_nbr,
                                                         mini, score_nodp, true, score_odpd,
                                                         mutex_miscellaneous);
                        } else {
                                score = score_nodp;
                                trace_pattern_and_scores(mram_info->nbr_len, current_read_nbr, ref_nbr,
                                                         mini, score_nodp, false, -1,
                                                         mutex_miscellaneous);
                        }

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
#ifdef DEBUG_RESULTS
                                        printf("[%u] nr=%u ix=%u num=%u ", me(), debug_nr, idx, request.num);
                                        printf("offset=%u seed=%u seq=%u score=%u\n",
                                               request.offset,
                                               ((uint32_t *) cached_coords_and_nbr)[0],
                                               ((uint32_t *) cached_coords_and_nbr)[1],
                                               score);
#endif /* DEBUG_RESULTS */
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

        stat_print(&stats);
        print_mram_stats(request_pool.stats_load, result_pool.stats_write, mutex_miscellaneous);

        mbox_send(&tasklet_stats, sizeof(tasklet_stats));
        return 0;
}

/**
 * @brief Common main: gather the seed mapping for this DPU, then start processing the requests on every tasklet.
 * @return 0 if everything worked fine
 */
int main() {
        barrier_t barrier = barrier_get(0);
        if (me() == 0) {
                load_mram_info();
                print_mram_info();
                request_pool_init(mram_info->io_offs);
                result_pool_init(MRAM_SIZE - MAX_DPU_RESULTS * sizeof(result_out_t));
#ifdef OPT_S
                if (((nb_bytes_to_syms(mram_info->nbr_len) + 1) * 24) >= 4096) {
                        printf("cannot run optimized code: symbol length is larger than mulub operation\n");
                        halt();
                }
#endif /* OPT_S */
                odpd_init(NB_RUNNING_TASKLETS, nb_bytes_to_syms(mram_info->nbr_len));
        }
        barrier_wait(barrier);

        /* For debugging purpose, one may reduce the number of operating tasklets. */
        if (me() < NB_RUNNING_TASKLETS) {
                return run_align();
        } else {
                return 0;
        }
}
