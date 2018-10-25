/*
 * Copyright (c) 2014-2018 - uPmem
 */

#include <stdlib.h>
#include <mdpu.h>
#include <stddef.h>
#include <upvc.h>
#include <string.h>
#include <dpus.h>

#include "dpus.h"
#include "dpulog.h"
#include "dispatch.h"
#include "parse_args.h"

static unsigned int nr_dpus_per_run;
unsigned int get_nr_dpus_per_run() { return nr_dpus_per_run;}

//#define LOG_DPUS
#ifdef LOG_DPUS
static dpu_logging_config_t logging_config = {
                                              .source = KTRACE,
                                              .destination_directory_name = "."
};
#define p_logging_config &logging_config
#define log_dpu(dpu, out) dpulog_read_for_dpu(dpu, out)
#else
#define p_logging_config NULL
#define log_dpu(dpu, log) do {} while(0)
#endif /* LOG_DPUS */

static struct dpu_param param = {
                                 .type = FUNCTIONAL_SIMULATOR,
                                 .profile = "",
                                 .on_boot = NULL,
                                 .logging_config = p_logging_config
};

void setup_dpus_for_target_type(target_type_t target_type) {
        switch (target_type) {
        case target_type_fpga:
                param.type = HW;
                param.profile = "chipId=6&enableDbgLib=true";
                nr_dpus_per_run = 8;
                break;

        default:
                param.type = FUNCTIONAL_SIMULATOR;
                nr_dpus_per_run = 1;
        }
}

devices_t dpu_try_alloc_for(unsigned int nr_dpus, const char *opt_program) {
        devices_t devices = (devices_t) malloc(sizeof(struct devices));
        if (devices == NULL) {
                ERROR_EXIT("*** could not allocate structure for devices - aborting");
        }
        devices->nr_dpus = nr_dpus;
        if (dpu_get_nr_of_dpus_for(&param, &(devices->nr_dpus_per_rank)) != DPU_API_SUCCESS) {
                ERROR_EXIT("*** could not guess the number of DPUs per rank - aborting");
        }
        devices->dpus = (dpu_t *) calloc(nr_dpus, sizeof(dpu_t));
        if (devices->dpus == NULL) {
                ERROR_EXIT("*** could not allocate array of DPUs - aborting");
        }
        unsigned int nr_ranks = (nr_dpus + devices->nr_dpus_per_rank - 1) / devices->nr_dpus_per_rank;
        devices->nr_ranks = nr_ranks;
        devices->ranks = (dpu_rank_t *) calloc(nr_ranks, sizeof(dpu_rank_t));
        if (devices->ranks == NULL) {
                ERROR_EXIT("*** could not allocate array of DPU ranks - aborting");
        }

        unsigned int each_rank, each_member, each_dpu;
        for (each_rank = each_dpu = 0; each_dpu < nr_dpus; each_rank++) {
                if (dpu_alloc(&param, &(devices->ranks[each_rank])) != DPU_API_SUCCESS) {
                        ERROR_EXIT("*** could not allocate rank number %u - aborting", each_rank);
                }
                for (each_member = 0;
                     (each_member < devices->nr_dpus_per_rank) && (each_dpu < nr_dpus); each_member++, each_dpu++) {
                        devices->dpus[each_dpu] = dpu_get_id(devices->ranks[each_rank], each_member);
                }
        }

        for (each_rank = 0; each_rank < devices->nr_ranks; each_rank++) {
                if (dpu_load_all(devices->ranks[each_rank], opt_program) != DPU_API_SUCCESS) {
                        ERROR_EXIT("*** could not load program on rank number %u - aborting", each_rank);
                }
        }

        return devices;
}

void dpu_try_load_mram_number(unsigned int mram_number, dpu_t dpu_id, devices_t devices, reads_info_t *reads_info) {
        mram_info_t *mram = mram_create(reads_info);
        mram_load(mram, mram_number);
        dpu_copy_to_individual(devices->dpus[dpu_id], (uint8_t *) mram, 0, mram->usage);
        mram_free(mram);
}

void dpu_try_free(devices_t devices) {
        if (devices->ranks != NULL) {
                unsigned int each_rank;
                for (each_rank = 0; each_rank < devices->nr_ranks; each_rank++) {
                        dpu_free(devices->ranks[each_rank]);
                }
                free(devices->ranks);
                devices->ranks = NULL;
        }
        if (devices->dpus != NULL) {
                free(devices->dpus);
                devices->dpus = NULL;
        }

        devices->nr_ranks = devices->nr_dpus = devices->nr_dpus_per_rank = 0;
}

uint64_t dpu_try_run(unsigned int dpu_id, devices_t devices) {
        // TODO: measure time
        uint64_t t0 = 0L; //dpu_time(dpu);
        if (dpu_boot_individual(devices->dpus[dpu_id], ASYNCHRONOUS) != DPU_API_SUCCESS) {
                ERROR("*** run failed on DPU number %u - aborting!", dpu_id);
                log_dpu(devices->dpus[dpu_id], stdout);
                exit(-8);
        }

        return t0;
}

bool dpu_try_check_status(unsigned int dpu_id, devices_t devices) {
        dpu_run_status_t run_status;
        if (dpu_get_individual_status(devices->dpus[dpu_id], &run_status) != DPU_API_SUCCESS) {
                ERROR_EXIT("*** could not get status from DPU number %u - aborting", dpu_id);
        }

        switch (run_status) {
        case DPU_STATUS_IDLE:
                return true;

        case DPU_STATUS_RUNNING:
                return false;

        case DPU_STATUS_ERROR:
                ERROR("*** DPU %u reported an error - aborting", dpu_id);
                log_dpu(devices->dpus[dpu_id], stdout);
                exit(-9);

        default:
                ERROR("*** could not get DPU status %u - aborting", dpu_id);
                log_dpu(devices->dpus[dpu_id], stdout);
                exit(-8);
        }
}

void dpu_try_write_dispatch_into_mram(unsigned int dpu_id,
                                      devices_t devices,
                                      unsigned int nr_reads,
                                      int8_t *reads,
                                      reads_info_t *reads_info) {
        mram_info_t mram;
        // Prepend the list of reads with information required by the DPU to properly align (see request_pool_init).
        struct {
                uint32_t nr_reads;
                uint32_t magic;
        } io_header = {
                       .nr_reads = nr_reads,
                       .magic = 0xcdefabcd
        };

        dpu_copy_from_individual(devices->dpus[dpu_id], (mram_addr_t) 0, (uint8_t *) &mram, sizeof(mram));
        unsigned int io_len = nr_reads * dispatch_read_len(reads_info);

        if (mram.io_offs + sizeof(io_header) + io_len >= MRAM_INPUT_SIZE) {
                ERROR_EXIT("*** will exceed MRAM limit if writing reads on DPU number %u - aborting!", dpu_id);
        }

        dpu_copy_to_individual(devices->dpus[dpu_id], (uint8_t *) &io_header, (mram_addr_t) (mram.io_offs),
                               sizeof(io_header));
        dpu_copy_to_individual(devices->dpus[dpu_id], (uint8_t *) reads, (mram_addr_t) (mram.io_offs + sizeof(io_header)), io_len);
        // Clear the output area
        uint8_t *reset_buffer = malloc(RESULT_AREA_LEN);
        memset(reset_buffer, -1, RESULT_AREA_LEN);
        dpu_copy_to_individual(devices->dpus[dpu_id], reset_buffer, (mram_addr_t) MRAM_INPUT_SIZE, RESULT_AREA_LEN);
        free(reset_buffer);
}

void dpu_try_log(unsigned int dpu_id, devices_t devices, uint64_t t0) {
        // TODO: get time
        uint64_t t1 = 0L; // dpu_time(dpu);

        double delta_in_secs = ((double) (t1 - t0)) / 750000000;
        printf("LOG DPU=%u TIME=%ld SEC=%.3f\n", dpu_id, t1 - t0, delta_in_secs);
        fflush(stdout);

        // Collect stats
        unsigned int each_tasklet;
        dpu_tasklet_stats_t tasklet_stats;
        for (each_tasklet = 0; each_tasklet < NR_TASKLET_PER_DPU; each_tasklet++) {
                if (dpu_tasklet_receive_individual(devices->dpus[dpu_id], (sysname_t) each_tasklet, 0,
                                                   (uint32_t *) &tasklet_stats,
                                                   sizeof(tasklet_stats)) != DPU_API_SUCCESS) {
                        ERROR_EXIT("*** could not read mailbox for DPU number %u tasklet %u", dpu_id, each_tasklet);
                }
                // TODO: show logs
                printf("LOG DPU=%u TID=%u REQ=%u\n", dpu_id, each_tasklet, tasklet_stats.nr_reqs);
                printf("LOG DPU=%u TID=%u NODP=%u\n", dpu_id, each_tasklet, tasklet_stats.nr_nodp_calls);
                printf("LOG DPU=%u TID=%u ODPD=%u\n", dpu_id, each_tasklet, tasklet_stats.nr_odpd_calls);
                printf("LOG DPU=%u TID=%u RESULTS=%u\n", dpu_id, each_tasklet, tasklet_stats.nr_results);
                printf("LOG DPU=%u TID=%u DATA_IN=%u\n", dpu_id, each_tasklet, tasklet_stats.mram_data_load);
                printf("LOG DPU=%u TID=%u RESULT_OUT=%u\n", dpu_id, each_tasklet, tasklet_stats.mram_result_store);
                printf("LOG DPU=%u TID=%u LOAD=%u\n", dpu_id, each_tasklet, tasklet_stats.mram_load);
                printf("LOG DPU=%u TID=%u STORE=%u\n", dpu_id, each_tasklet, tasklet_stats.mram_store);
        }
        fflush(stdout);

        log_dpu(devices->dpus[dpu_id], stdout);
}

dpu_result_out_t *dpu_try_get_results(unsigned int dpu_id, devices_t devices) {
        dpu_result_out_t *result_buffer = (dpu_result_out_t *) malloc(RESULT_AREA_LEN);
        dpu_copy_from_individual(devices->dpus[dpu_id], (mram_addr_t) MRAM_INPUT_SIZE, (uint8_t *) result_buffer,
                                 RESULT_AREA_LEN);
        return result_buffer;
}

void dpu_try_backup_mram(unsigned int tid, devices_t devices, const char *file_name) {
        printf("saving DPU %u MRAM into '%s'\n", tid, file_name);
        FILE *f = fopen(file_name, "wb");
        if (f != NULL) {
                uint8_t *mram = (uint8_t *) malloc(1 << 26);
                dpu_copy_from_individual(devices->dpus[tid], 0, mram, 1 << 26);
                fwrite(mram, 1, 1 << 26, f);
                free(mram);
                fclose(f);
        } else {
                (void) fprintf(stderr, "WARNING: failed to backup DPU %u MRAM into '%s'\n", tid, file_name);
        }
}

