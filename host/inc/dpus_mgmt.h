/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DPUS_H__
#define __INTEGRATION_DPUS_H__

#include <dpu.h>
#include <dpulog.h>

#include "parse_args.h"
#include "upvc.h"
#include "dispatch.h"

#include "common.h"

/**
 * @brief General structure describing the DPUs involved in a run.
 */
typedef struct {
        unsigned int nb_dpus_per_rank;
        unsigned int nb_ranks_per_run;
        unsigned int nb_ranks;
        dpu_rank_t *ranks;
        unsigned int nb_dpus;
        dpu_t *dpus;
        mram_info_t *mram_info;
} devices_t;

/**
 * @brief Preliminary setup operation that defines the profile of allocated DPUs according to the target type.
 *
 * @param target_type  The type of targeted DPUs (hsim or fpga).
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
 * @param dpu_id       The DPU number.
 * @param devices      Available devices.
 * @param mram         The mram content.
 */
void dpu_try_write_mram(dpu_t dpu_id,
                        devices_t *devices,
                        mram_info_t *mram);

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
 * @brief Checks the status of a rank of DPUs.
 *
 * @param rank_id  The rank number.
 * @param devices  Available devices.
 *
 * @return Whether the DPU has finished its run.
 */
bool dpu_try_check_status(unsigned int rank_id, devices_t *devices);

/**
 * @brief Stores the dispatch parameters into a DPU, raises an exception if the MRAM limit is reached.
 *
 * @param rank_id     The rank number.
 * @param dpu_offset  The offset in DPUs to find the mram.
 * @param devices     List of available devices.
 * @param dispatch    The structure containing the dispatching of the reads into the DPUs.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
void dpu_try_write_dispatch_into_mram(unsigned int rank_id,
                                      unsigned int dpu_offset,
                                      devices_t *devices,
                                      dispatch_request_t *dispatch,
                                      reads_info_t *reads_info);

/**
 * @brief Reads the result area of a DPU, raises an exception if something went wrong with the DPU.
 *
 * @param dpu_id   The DPU number.
 * @param devices  List of available devices.
 * @param result_buffer The table of DPU results to store the results (last one having a request number equal to -1).
 */
void dpu_try_get_results(unsigned int dpu_id, devices_t *devices, dpu_result_out_t *result_buffer);

/**
 * @brief Reports the log of a DPU and the execution time.
 *
 * @param dpu_id   The DPU number.
 * @param devices  Available devices.
 */
void dpu_try_log(unsigned int dpu_id, devices_t *devices);

/**
 * @brief Makes a snapshot of an MRAM into a given file.
 *
 * @param dpu_id     The DPU number.
 * @param devices    List of available devices.
 * @param file_name  Where to save data.
 */
void dpu_try_backup_mram(unsigned int tid, devices_t *devices, const char *file_name);

#endif /* __INTEGRATION_DPUS_H__ */
