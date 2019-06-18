/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <assert.h>
#include <stdint.h>

/* #define STATS_ON */

#define MRAM_SIZE (64 << 20)

#define NB_TASKLET_PER_DPU_LOG2 4
#define NB_TASKLET_PER_DPU (1 << NB_TASKLET_PER_DPU_LOG2)

#define ALIGN_DPU(val) (((val) + 7) & ~7)

#define MAX_DPU_REQUEST (1 << 16)
#define MAX_DPU_RESULTS (1 << 16)
#define MAX_RESULTS_PER_READ (1 << 10)

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
#define MRAM_INFO_ADDR(mram_offset) (ALIGN_DPU(mram_offset))
#define MRAM_INFO_READ(mram_info, mram_offset)                                                                                   \
    do {                                                                                                                         \
        mram_read16(MRAM_INFO_ADDR(mram_offset), mram_info);                                                                     \
    } while (0)
_Static_assert(sizeof(mram_info_t) == 16, "mram_info_t size changed (make sure that MRAM_INFO_READ changed as well)");

/**
 * @brief Coordonates of the read that matched in the reference genome.
 */
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
 *
 * @var num  Number that associate an input read with a request.
 * @var score Best score of the read with a read of the reference genome.
 * @var coord Coordinate of the read that matched in the reference genome.
 */
typedef struct {
    int32_t num;
    uint32_t score;
    dpu_result_coord_t coord;
} dpu_result_out_t;
#define DPU_RESULT_VAR m_dpu_result

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
    uint64_t nodp_time;
    uint64_t odpd_time;
} dpu_tasklet_stats_t;
#define DPU_TASKET_STATS_VAR m_dpu_tasklet_stats

typedef uint64_t dpu_compute_time_t;
#define DPU_COMPUTE_TIME_VAR m_dpu_compute_time

#define DPU_INPUTS_ADDR(mram_offset) (ALIGN_DPU(MRAM_INFO_ADDR(mram_offset) + sizeof(mram_info_t)))
#define DPU_INPUTS_SIZE(mram_offset) (MRAM_SIZE - DPU_INPUTS_ADDR(mram_offset))

/**
 * @brief Information on the requests reads to a DPU.
 */
typedef struct {
    uint32_t nb_reads;
    uint32_t magic;
} request_info_t;

#define DPU_REQUEST_INFO_ADDR(mram_info, mram_offset) (ALIGN_DPU(DPU_INPUTS_ADDR(mram_offset) + (mram_info)->total_nbr_size))

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

#define DPU_REQUEST_ADDR(mram_info, mram_offset)                                                                                 \
    (ALIGN_DPU(DPU_REQUEST_INFO_ADDR(mram_info, mram_offset) + sizeof(request_info_t)))
#define DPU_REQUEST_SIZE(size_nbr) (ALIGN_DPU(sizeof(dpu_request_t) + (size_nbr)))

#endif /* __COMMON_H__ */
