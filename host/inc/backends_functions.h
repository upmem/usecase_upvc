#ifndef __BACKENDS_FUNCTIONS_H__
#define __BACKENDS_FUNCTIONS_H__

#include <semaphore.h>

/**
 * @brief Specific function whether the application run on simulation or on DPU (hsim of fpga).
 */
typedef struct backends_functions_struct {
        void (*run_dpu)(dispatch_request_t *,
                        devices_t *,
                        unsigned int,
                        unsigned int,
                        unsigned int,
                        sem_t *,
                        sem_t *,
                        times_ctx_t *,
                        reads_info_t *);
        void (*add_seed_to_requests)(dispatch_request_t *, int, int, index_seed_t *, int8_t *, reads_info_t *);

        vmi_t *(*init_vmis)(unsigned int);
        void (*free_vmis)(vmi_t *, unsigned int, unsigned int *, reads_info_t *);
        void (*write_vmi)(vmi_t *, unsigned int, unsigned int, int8_t *, uint64_t, reads_info_t *);

        void (*init_backend)(unsigned int *,
                             devices_t **,
                             unsigned int,
                             const char *,
                             index_seed_t ***,
                             unsigned int,
                             genome_t *,
                             reads_info_t *,
                             times_ctx_t *,
                             struct backends_functions_struct *);
        void (*free_backend)(devices_t *, unsigned int);
        void (*load_mram)(unsigned int, unsigned int, devices_t *, reads_info_t *, times_ctx_t *);
} backends_functions_t;

#endif /* __BACKENDS_FUNCTIONS_H__ */
