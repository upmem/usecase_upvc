#ifndef __DPU_BACKEND_H__
#define __DPU_BACKEND_H__

#include "genome.h"
#include "vmi.h"
#include "upvc.h"
#include "index.h"
#include "dispatch.h"
#include "backends_functions.h"

/**
 * @brief get index of seed from files on disk
 */
index_seed_t **get_index_seed_dpu(unsigned int nb_dpu,
                                  __attribute__((unused)) genome_t *ref_genome,
                                  reads_info_t *reads_info,
                                  __attribute__((unused)) times_ctx_t *times_ctx,
                                  __attribute__((unused)) backends_functions_t *backends_funcions);

/**
 * @brief Init the Virtual Memory Images structure needed to produce the mram images.
 */
vmi_t *init_vmis_dpu(unsigned int nb_dpu);

/**
 * @brief Free the Virtual Memory Images.
 */
void free_vmis_dpu(vmi_t *vmis, unsigned int nb_dpu, unsigned int *nb_neighbours, reads_info_t *reads_info);

/**
 * @brief Write data in a Virtual Memory Image.
 */
void write_vmi_dpu(vmi_t *vmis, unsigned int dpuno, unsigned int k, int8_t *nbr, uint64_t coords, reads_info_t *reads_info);

/**
 * @brief Add a seed to a request.
 */
void add_seed_to_dpu_requests(dispatch_request_t *requests,
                              int num_dpu,
                              int num_read,
                              __attribute__((unused)) int nb_read_written,
                              index_seed_t *seed,
                              int8_t *nbr,
                              reads_info_t *reads_info);

/**
 * @brief Compute one pass on DPUs.
 */
void run_on_dpu(dispatch_t dispatch,
                devices_t *devices,
                unsigned int nb_dpu,
                times_ctx_t *times_ctx,
                reads_info_t *reads_info);

#endif /* __DPU_BACKEND_H__ */
