#ifndef __SIMU_BACKEND_H__
#define __SIMU_BACKEND_H__

#include <semaphore.h>

#include "genome.h"
#include "vmi.h"
#include "upvc.h"
#include "index.h"
#include "dispatch.h"
#include "backends_functions.h"

/**
 * @brief Does nothing for simulation
 */
vmi_t *init_vmis_simulation( unsigned int nb_dpu);

/**
 * @brief Does nothing for simulation
 */
void free_vmis_simulation( vmi_t *vmis,
                           unsigned int nb_dpu,
                           unsigned int *nb_neighbours,
                           reads_info_t *reads_info);

/**
 * @brief Write neighbour and coordinate in structure representing the DPUs in simulation mode
 */
void write_vmi_simulation( vmi_t *vmis,
                          unsigned int dpuno,
                          unsigned int align_idx,
                          int8_t *nbr,
                          dpu_result_coord_t coord,
                          reads_info_t *reads_info);

/**
 * @brief Read the reference genome and create the index seed with it.
 */
index_seed_t **get_index_seed_simulation();

/**
 * @brief Write seed information in the structure representing the DPU list of requests.
 */
void add_seed_to_simulation_requests(dispatch_request_t *requests,
                                     int num_read,
                                     int nb_read_written,
                                     index_seed_t *seed,
                                     int8_t *nbr,
                                     reads_info_t *reads_info);

/**
 * @brief Compute one pass in simulation mode.
 */
void run_dpu_simulation(dispatch_request_t *dispatch,
                        devices_t *devices,
                        unsigned int dpu_offset,
                        unsigned int rank_id,
                        unsigned int nb_pass,
                        sem_t *dispatch_free_sem,
                        sem_t *acc_wait_sem,
                        times_ctx_t *times_ctx,
                        reads_info_t *reads_info);

/**
 * @brief Index the reference genome.
 */
void init_backend_simulation(unsigned int *nb_rank,
                             devices_t **devices,
                             unsigned int nb_dpu_per_run,
                             const char *dpu_binary,
                             index_seed_t ***index_seed,
                             unsigned int nb_dpu,
                             genome_t *ref_genome,
                             reads_info_t *reads_info,
                             times_ctx_t *times_ctx,
                             backends_functions_t *backends_functions);

/**
 * @brief Free structure used by simulation.
 */
void free_backend_simulation(devices_t *devices, unsigned int nb_dpu);

/**
 * @brief Does nothing in sumation mode.
 */
void load_mram_simulation(unsigned int dpu_offset,
                          unsigned int rank_id,
                          devices_t *devices,
                          reads_info_t *reads_info,
                          times_ctx_t *times_ctx);

#endif /* __SIMU_BACKEND_H__ */
