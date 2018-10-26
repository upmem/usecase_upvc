/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DPUS_H__
#define __INTEGRATION_DPUS_H__

#include <dpu.h>
#include <dpulog.h>
#include "parse_args.h"
#include "upvc.h"

/**
 * @brief Bumber of tasklets handled by each DPU.
 */
#define NB_TASKLET_PER_DPU_LOG2 4
#define NB_TASKLET_PER_DPU (1 << NB_TASKLET_PER_DPU_LOG2)

/**
 * @brief General structure describing the DPUs involved in a run.
 */
typedef struct devices {
        unsigned int nb_dpus_per_rank;
        unsigned int nb_ranks;
        dpu_rank_t *ranks;
        unsigned int nb_dpus;
        dpu_t *dpus;
} *devices_t;

/**
 * @brief Results returned by DPUs.
 */
typedef struct {
        int num;
        uint32_t seq_nr;
        uint32_t seed_nr;
        unsigned int score;
} dpu_result_out_t;

/**
 * @brief Statistics reported by individual tasklets.
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
 * @brief Get the number of DPUs to be run per run.
 *
 * @return The number of DPUs to be run per run.
 */
unsigned int get_nb_dpus_per_run();

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
devices_t dpu_try_alloc_for(unsigned int nb_dpus, const char *opt_program);

/**
 * @brief Initializes the contents of an MRAM.
 *
 * @param mram_number  It identifies the MRAM to load into a DPU.
 * @param dpu_id       The DPU number.
 * @param devices      Available devices.
 * @param reads_info   Information on the size of the seed and the neighbour.
 */
void dpu_try_load_mram_number(unsigned int mram_number, dpu_t dpu_id, devices_t devices, reads_info_t *reads_info);

/**
 * @brief Frees the allocated DPUs.
 *
 * @param devices  The list of allocated devices.
 */
void dpu_try_free(devices_t devices);

/**
 * @brief Boots a DPU, raises an exception if the execution fails.
 *
 * The DPU is simply booted. Further calls to dpu_try_get_status will tell whether the execution
 * is complete.
 *
 * @param dpu_id   The DPU number.
 * @param devices  Available devices.
 *
 * @return DPU time, before booting.
 */
uint64_t dpu_try_run(unsigned int dpu_id, devices_t devices);

/**
 * @brief Checks the status of a DPU.
 *
 * @param dpu_id   The DPU number.
 * @param devices  Available devices.
 *
 * @return Whether the DPU has finished its run.
 */
bool dpu_try_check_status(unsigned int dpu_id, devices_t devices);

/**
 * @brief Stores the dispatch parameters into a DPU, raises an exception if the MRAM limit is reached.
 *
 * @param dpu_id      The DPU number.
 * @param devices     List of available devices.
 * @param nb_reads    How many reads to write.
 * @param reads       A table of nb_reads requests.
 * @param reads_info  Information on the size of the seed and the neighbour.
 */
void dpu_try_write_dispatch_into_mram(unsigned int dpu_id,
                                      devices_t devices,
                                      unsigned int nb_reads,
                                      int8_t *reads,
                                      reads_info_t *reads_info);

/**
 * @brief Reads the result area of a DPU, raises an exception if something went wrong with the DPU.
 *
 * @param dpu_id   The DPU number.
 * @param devices  List of available devices.
 *
 * @return A table of DPU results, the last one having a request number equal to -1, which must be freed when done.
 */
dpu_result_out_t *dpu_try_get_results(unsigned int dpu_id, devices_t devices);

/**
 * @brief Reports the log of a DPU and the execution time.
 *
 * @param dpu_id   The DPU number.
 * @param devices  Available devices.
 * @param t0       The initial DPU time, reported by dpu_try_run.
 */
void dpu_try_log(unsigned int dpu_id, devices_t devices, uint64_t t0);

/**
 * @brief Makes a snapshot of an MRAM into a given file.
 *
 * @param dpu_id     The DPU number.
 * @param devices    List of available devices.
 * @param file_name  Where to save data.
 */
void dpu_try_backup_mram(unsigned int tid, devices_t devices, const char *file_name);

#endif /* __INTEGRATION_DPUS_H__ */
