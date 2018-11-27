#ifndef __SIMU_BACKEND_H__
#define __SIMU_BACKEND_H__

#include "genome.h"
#include "vmi.h"
#include "upvc.h"
#include "index.h"
#include "dispatch.h"
#include "backends_functions.h"

/**
 * @brief Does nothing for simulation
 */
vmi_t *init_vmis_simulation(__attribute__((unused)) unsigned int nb_dpu);

/**
 * @brief Does nothing for simulation
 */
void free_vmis_simulation(__attribute__((unused)) vmi_t *vmis,
                          __attribute__((unused)) unsigned int nb_dpu,
                          __attribute__((unused)) unsigned int *nb_neighbours,
                          __attribute__((unused)) reads_info_t *reads_info);

/**
 * @brief Write neighbour and coordinate in structure representing the DPUs in simulation mode
 */
void write_vmi_simulation(__attribute__((unused)) vmi_t *vmis,
                          unsigned int dpuno,
                          unsigned int align_idx,
                          int8_t *nbr,
                          long coords,
                          reads_info_t *reads_info);

/**
 * @brief Read the reference genome and create the index seed with it.
 */
index_seed_t **get_index_seed_simulation(unsigned int nb_dpu,
                                         genome_t *ref_genome,
                                         reads_info_t *reads_info,
                                         times_ctx_t *times_ctx,
                                         backends_functions_t *backends_functions);

/**
 * @brief Write seed information in the structure representing the DPU list of requests.
 */
void add_seed_to_simulation_requests(__attribute__((unused)) dispatch_request_t *requests,
                                     __attribute__((unused)) int num_dpu,
                                     int num_read,
                                     int nb_read_written,
                                     index_seed_t *seed,
                                     int8_t *nbr,
                                     reads_info_t *reads_info);

/**
 * @brief Compute one pass in simulation mode
 */
void run_dpu_simulation(__attribute__((unused)) dispatch_t dispatch,
                        __attribute__((unused)) devices_t devices,
                        unsigned int nb_dpu,
                        times_ctx_t *times_ctx,
                        reads_info_t *reads_info);

/**
 * @brief Does nothing in simulation
 */
devices_t init_devices_simulation(unsigned int nb_dpu,
                                  const char *dpu_binary);

/**
 * @brief Does nothing in simulation
 */
void free_devices_simulation(devices_t devices);

#endif /* __SIMU_BACKEND_H__ */
