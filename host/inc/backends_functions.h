#ifndef __BACKENDS_FUNCTIONS_H__
#define __BACKENDS_FUNCTIONS_H__

/**
 * @brief Specific function whether the application run on simulation or on DPU (hsim of fpga).
 */
typedef struct backends_functions_struct {
        void (*run_dpu)(dispatch_t, devices_t *, unsigned int, times_ctx_t *, reads_info_t *);
        void (*add_seed_to_requests)(dispatch_request_t *, int, int, int, index_seed_t *, int8_t *, reads_info_t *);
        index_seed_t **(*get_index_seed)(unsigned int,
                                         genome_t *,
                                         reads_info_t *,
                                         times_ctx_t *,
                                         struct backends_functions_struct *);

        vmi_t *(*init_vmis)(unsigned int);
        void (*free_vmis)(vmi_t *, unsigned int, unsigned int *, reads_info_t *);
        void (*write_vmi)(vmi_t *, unsigned int, unsigned int, int8_t *, uint64_t, reads_info_t *);

        devices_t *(*init_devices)(unsigned int, const char *);
        void (*free_devices)(devices_t *);
} backends_functions_t;

#endif /* __BACKENDS_FUNCTIONS_H__ */
