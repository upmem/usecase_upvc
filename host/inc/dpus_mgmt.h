/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DPUS_H__
#define __INTEGRATION_DPUS_H__

#include <dpu.h>
#include <dpu_log.h>

#include "dispatch.h"
#include "parse_args.h"
#include "upvc.h"

#include "common.h"

/**
 * @brief General structure describing the DPUs involved in a run.
 */
typedef struct {
    unsigned int nb_ranks_per_run;
    unsigned int *nb_dpus_per_rank;
    unsigned int *rank_mram_offset;
    struct dpu_rank_t **ranks;
    unsigned int nb_dpus;
    pthread_mutex_t log_mutex;
    FILE *log_file;
    struct dpu_symbol_t mram_info, mram_request_info, mram_requests, mram_compute_time, mram_result, mram_tasklet_stats,
        mram_available, mram_results_checksum;
} devices_t;

/**
 * @brief Preliminary setup operation that defines the profile of allocated DPUs according to the target type.
 *
 * @param target_type  The type of targeted DPUs (fsim or fpga).
 */
void setup_dpus_for_target_type(target_type_t target_type);

/**
 * @brief Creates the devices to handle a given number of DPUs.
 *
 * DPUs are created, loaded with the same program.
 *
 * @param nb_dpus      How many DPUs are required.
 * @param opt_program  Program loaded into every DPU.
 *
 * @return The allocated devices structure, or NULL in case of error.
 */
devices_t *dpu_try_alloc_for(unsigned int nb_dpus, const char *opt_program);

/**
 * @brief Initializes the contents of an MRAM.
 *
 * @param rank_id      The rank number.
 * @param devices      Available devices.
 * @param mram         The mram content.
 * @param mram_size    The size of the mrams.
 * @param delta_neighbour   The delta to apply to the UPVC DPU algorithm.
 */
void dpu_try_write_mram(unsigned int rank_id, devices_t *devices, uint8_t **mram, size_t *mram_size, int delta_neighbour);

/**
 * @brief Frees the allocated DPUs.
 *
 * @param devices  The list of allocated devices.
 */
void dpu_try_free(devices_t *devices);

/**
 * @brief Boots a rank of DPUs, raises an exception if the execution fails.
 *
 * The DPUs are simply booted. Further calls to dpu_try_get_status will tell whether the execution
 * is complete.
 *
 * @param rank_id  The rank number.
 * @param devices  Available devices.
 */
void dpu_try_run(unsigned int rank_id, devices_t *devices);

/**
 * @brief Stores the dispatch parameters into a DPU, raises an exception if the MRAM limit is reached.
 *
 * @param rank_id     The rank number.
 * @param dpu_offset  The offset in DPUs to find the mram.
 * @param devices     List of available devices.
 * @param dispatch    The structure containing the dispatching of the reads into the DPUs.
 */
void dpu_try_write_dispatch_into_mram(
    unsigned int rank_id, unsigned int dpu_offset, devices_t *devices, dispatch_request_t *dispatch, uint64_t *requests_checksum);

/**
 * @brief Reads the result area of a DPU, raises an exception if something went wrong with the DPU.
 *        Also print the log (statistic).
 *
 * @param rank_id  The rank number.
 * @param dpu_offset The offset in the DPUs.
 * @param devices  List of available devices.
 * @param result_buffer The table of DPU results to store the results (last one having a request number equal to -1).
 */
void dpu_try_get_results_and_log(unsigned int rank_id, unsigned int dpu_offset, devices_t *devices,
    dpu_result_out_t **result_buffer, uint64_t *results_checksum);

#endif /* __INTEGRATION_DPUS_H__ */
