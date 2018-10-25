/*
 * Copyright (c) 2014-2018 - uPmem
 */

#ifndef INTEGRATION_DPUS_H
#define INTEGRATION_DPUS_H

#include <dpu.h>
#include <dpulog.h>
#include "parse_args.h"
#include "upvc.h"

/**
 * @brief number of tasklets handled by each DPU
 */
#define NR_TASKLET_PER_DPU_LOG2 4
#define NR_TASKLET_PER_DPU (1 << NR_TASKLET_PER_DPU_LOG2)

unsigned int get_nr_dpus_per_run();

/**
 * @brief General structure describing the DPUs involved in a run.
 */
typedef struct devices {
    unsigned int nr_dpus_per_rank;
    unsigned int nr_ranks;
    dpu_rank_t *ranks;
    unsigned int nr_dpus;
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
    uint32_t nr_reqs;
    uint32_t nr_nodp_calls;
    uint32_t nr_odpd_calls;
    uint32_t nr_results;
    uint32_t mram_data_load;
    uint32_t mram_result_store;
    uint32_t mram_load;
    uint32_t mram_store;
} dpu_tasklet_stats_t;

/**
 * @brief Preliminary setup operation that defines the profile of allocated DPUs according to the target type
 * @param target_type the type of targeted DPUs
 */
void setup_dpus_for_target_type(target_type_t target_type);

/**
 * @brief Creates the devices to handle a given number of DPUs
 * @param nr_dpus how many DPUs are required
 * @param opt_program program loaded into every DPU (NULL for the default one)
 * @return The allocated devices structure, or NULL in case of error
 *
 * DPUs are created, loaded with the same program.
 */
devices_t dpu_try_alloc_for(unsigned int nr_dpus, const char *opt_program);

/**
 * @brief Initializes the contents of an MRAM
 * @param mram_number identifies the MRAM to load into a DPU
 * @param dpu_id the DPU number
 * @param devices available devices
 * @param reads_info are information on the seed and neighbour size
 */
void dpu_try_load_mram_number(unsigned int mram_number, dpu_t dpu_id, devices_t devices, reads_info_t *reads_info);

/**
 * @brief Frees the allocated DPUs
 * @param devices the list of allocated devices
 */
void dpu_try_free(devices_t devices);

/**
 * @brief Boots a DPU, raises an exception if the execution fails.
 *
 * The DPU is simply booted. Further calls to dpu_try_get_status will tell whether the execution
 * is complete.
 *
 * @param dpu_id the DPU number
 * @param devices available devices
 * @return DPU time, before booting
 */
uint64_t dpu_try_run(unsigned int dpu_id, devices_t devices);

/**
 * @brief Checks the status of a DPU
 * @param dpu_id the DPU number
 * @param devices available devices
 * @return Whether the DPU has finished its run
 */
bool dpu_try_check_status(unsigned int dpu_id, devices_t devices);

/**
 * @brief Stores the dispatch parameters into a DPU, raises an exception if the MRAM limit is reached.
 * @param dpu_id   the DPU number
 * @param devices  list of available devices
 * @param nr_reads how many reads to write
 * @param reads    a table of nr_reads requests
 * @param reads_info are information on the seed and neighbour size
 */
void dpu_try_write_dispatch_into_mram(unsigned int dpu_id,
                                      devices_t devices,
                                      unsigned int nr_reads,
                                      int8_t *reads,
                                      reads_info_t *reads_info);

/**
 * @brief Reads the result area of a DPU, raises an exception if something went wrong with the DPU.
 * @param dpu_id   the DPU number
 * @param devices  list of available devices
 * @return A table of DPU results, the last one having a request number equal to -1, which must be freed when done.
 */
dpu_result_out_t *dpu_try_get_results(unsigned int dpu_id, devices_t devices);

/**
 * @brief Reports the log of a DPU and the execution time
 * @param dpu_id   the DPU number
 * @param devices  available devices
 * @param t0       the initial DPU time, reported by dpu_try_run
 */
void dpu_try_log(unsigned int dpu_id, devices_t devices, uint64_t t0);

/**
 * @brief Makes a snapshot of an MRAM into a given file
 * @param dpu_id    the DPU number
 * @param devices   list of available devices
 * @param file_name where to save data
 */
void dpu_try_backup_mram(unsigned int tid, devices_t devices, const char *file_name);

#endif //INTEGRATION_DPUS_H
