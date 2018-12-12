#ifndef __COMMON_H__
#define __COMMON_H__

#include "assert.h"

#define STATS_ON

#define MRAM_SIZE (64 << 20)

#define NB_TASKLET_PER_DPU_LOG2 4
#define NB_TASKLET_PER_DPU (1 << NB_TASKLET_PER_DPU_LOG2)

#define ALIGN_DPU(val) (((val) + 7) & ~7)

#define MAX_DPU_REQUEST      (1 << 19)
#define MAX_DPU_RESULTS      (1 << 19)
#define MAX_RESULTS_PER_READ (1 << 10)

/**
 * DPU MRAM Memory Layout
 *
 * 0 MRAM_INFO_ADDR
 * | sizeof(mram_info_t)
 * |
 * | DPU_INPUTS_ADDR
 * | mram_info.total_nbr_size
 * |
 * | DPU_REQUEST_INFO_ADDR
 * | sizeof(request_info_t)
 * |
 * | DPU_REQUEST_ADDR
 * | MAX_DPU_REQUEST * DPU_REQUEST_SIZE
 * |
 * | EMPTY SPACE
 * |
 * | DPU_TASKLET_COMPUTE_TIME_ADDR
 * | DPU_TASKLET_COMPUTE_TIME_SIZE
 * |
 * | DPU_TASKLET_STATS_ADDR
 * | DPU_TASKLET_STATS_SIZE
 * |
 * | DPU_SWAP_RESULT_ADDR
 * | DPU_SWAP_RESULT_SIZE
 * |
 * | DPU_RESULT_ADDR
 * | DPU_RESULT_SIZE
 * MRAM_SIZE
 *
 */

/**
 * @brief A snapshot of the MRAM.
 *
 * @var total_nbr_size  Size of all the reads of the reference genome.
 * @var nb_nbr          Total number of neighbours stored in this MRAM.
 * @var nbr_len         Length of a neighbour, in bytes.
 * @var delta           Delta to apply to nodp and odpd comparison depending on the round.
 */
typedef struct {
        uint32_t total_nbr_size;
        uint32_t nb_nbr;
        uint32_t nbr_len;
        uint32_t delta;
} mram_info_t;

#define MRAM_INFO_ADDR ((mram_addr_t)0)
#define MRAM_INFO_READ(addr, mram_info) do { mram_read16(addr, mram_info); } while(0)
_Static_assert(sizeof(mram_info_t) == 16 ,"mram_info_t size changed (make sure that MRAM_INFO_READ changed as well)");

typedef struct {
        union {
                uint64_t coord;
                struct {
                        uint32_t seed_nr;
                        uint32_t seq_nr;
                };
        };
} dpu_result_coord_t;
/**
 * @brief One result produced for one read
 */
typedef struct {
        int32_t num;
        uint32_t score;
        dpu_result_coord_t coord;
} dpu_result_out_t;

#define DPU_RESULT_SIZE (MAX_DPU_RESULTS * sizeof(dpu_result_out_t))
#define DPU_RESULT_ADDR (MRAM_SIZE - DPU_RESULT_SIZE)
#define DPU_RESULT_WRITE(res, addr) do { mram_write16(res, addr); } while(0)
_Static_assert(sizeof(dpu_result_out_t) == 16, "dpu_result_out_t size changed (make sure that DPU_RESULT_WRITE changed as well)");

#define DPU_SWAP_RESULT_SIZE (NB_TASKLET_PER_DPU * MAX_RESULTS_PER_READ * sizeof(dpu_result_out_t))
#define DPU_SWAP_RESULT_ADDR (DPU_RESULT_ADDR - DPU_SWAP_RESULT_SIZE)

/**
 * @brief stats reported by every tasklet
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

#define DPU_TASKLET_STATS_SIZE (sizeof(dpu_tasklet_stats_t) * NB_TASKLET_PER_DPU)
#define DPU_TASKLET_STATS_ADDR (DPU_SWAP_RESULT_ADDR - DPU_TASKLET_STATS_SIZE)
#define DPU_TASKLET_STATS_WRITE(res, addr) do { mram_write32(res, addr); } while(0)
_Static_assert(sizeof(dpu_tasklet_stats_t) == 32,
               "dpu_tasklet_stats_t size changed (make sure that DPU_TASKLET_STATS_WRITE changed as well)");

typedef uint64_t dpu_compute_time_t;
#define DPU_COMPUTE_TIME_SIZE (sizeof(dpu_compute_time_t))
#define DPU_COMPUTE_TIME_ADDR (DPU_TASKLET_STATS_ADDR - DPU_COMPUTE_TIME_SIZE)
#define DPU_COMPUTE_TIME_WRITE(res, addr) do { mram_write8(res, addr); } while(0)
_Static_assert(sizeof(dpu_compute_time_t) == 8,
               "dpu_compute_time_t size changed (make sure that DPU_COMPUTE_TIME_WRITE changed as well)");

#define DPU_INPUTS_ADDR (ALIGN_DPU(MRAM_INFO_ADDR + sizeof(mram_info_t)))
#define DPU_INPUTS_SIZE (DPU_COMPUTE_TIME_ADDR - DPU_INPUTS_ADDR)

/**
 * @brief Information on the requests reads to a DPU.
 */
typedef struct {
        uint32_t nb_reads;
        uint32_t magic;
} request_info_t;

#define DPU_REQUEST_INFO_ADDR(mram_info) (ALIGN_DPU(DPU_INPUTS_ADDR + (mram_info)->total_nbr_size))

#define REQUEST_INFO_READ(addr, request_info) do { mram_read8(addr, request_info); } while(0)
_Static_assert(sizeof(request_info_t) == 8, "request_info_t size changed (make sure that REQUEST_INFO_READ changed as well)");

/**
 * @brief Structure representing one request to a DPU, resulting from the dispatching of reads.
 *
 * Such a request basically contains all the information for one read.
 *
 * @var offset  The 1st neighbour address.
 * @var count   The number of neighbours.
 * @var num     A reference number to the original request.
 */
typedef struct {
        uint32_t offset;
        uint32_t count;
        uint32_t num;
} dpu_request_t;

#define DPU_REQUEST_ADDR(mram_info) (ALIGN_DPU(ALIGN_DPU(DPU_INPUTS_ADDR + (mram_info)->total_nbr_size) + sizeof(request_info_t)))
#define DPU_REQUEST_SIZE(size_nbr) (ALIGN_DPU(sizeof(dpu_request_t) + (size_nbr)))

#endif /* __COMMON_H__ */
