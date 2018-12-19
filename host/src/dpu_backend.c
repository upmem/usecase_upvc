#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <semaphore.h>

#include "mram_dpu.h"
#include "dpus_mgmt.h"
#include "upvc_dpu.h"
#include "genome.h"
#include "index.h"
#include "dispatch.h"
#include "upvc.h"
#include "vmi.h"
#include "backends_functions.h"
#include "common.h"
#include "parse_args.h"

vmi_t *init_vmis_dpu(unsigned int nb_dpu)
{
        vmi_t *vmis = (vmi_t *) calloc(nb_dpu, sizeof(vmi_t));
        char vmi_name[8];
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                (void) sprintf(vmi_name, "%04u", dpuno);
                vmi_create(vmi_name, vmis + dpuno);
        }
        return vmis;
}

static void dump_mdpu_images_into_mram_files(vmi_t *vmis, unsigned int *nbr_counts, unsigned int nb_dpu, reads_info_t *reads_info)
{
        mram_info_t *mram_image = (mram_info_t *)malloc(MRAM_SIZE);
        printf("Creating MRAM images\n");
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
                vmi_t *this_vmi = vmis + each_dpu;
                mram_copy_vmi(mram_image, this_vmi, nbr_counts[each_dpu], reads_info);
                mram_save(mram_image, each_dpu);
        }
        free(mram_image);
}

void free_vmis_dpu(vmi_t *vmis, unsigned int nb_dpu, unsigned int *nb_neighbours, reads_info_t *reads_info)
{
        dump_mdpu_images_into_mram_files(vmis, nb_neighbours, nb_dpu, reads_info);
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                vmi_delete(vmis + dpuno);
        }
        free(vmis);
}

void write_vmi_dpu(vmi_t *vmis, unsigned int dpuno, unsigned int k, int8_t *nbr, uint64_t coords, reads_info_t *reads_info)
{
        unsigned int size_neighbour_in_bytes = reads_info->size_neighbour_in_bytes;
        unsigned int out_len = ALIGN_DPU(sizeof(uint64_t) + size_neighbour_in_bytes);
        uint64_t temp_buff[out_len / sizeof(uint64_t)];
        memset(temp_buff, 0, out_len);
        temp_buff[0] = coords;
        memcpy(&temp_buff[1], nbr, (size_t) size_neighbour_in_bytes);
        vmi_write(vmis + dpuno, k * out_len, temp_buff, out_len);
}

static void dispatch_request_add(dispatch_request_t *reads,
                                 unsigned int offset,
                                 unsigned int count,
                                 unsigned int num,
                                 int8_t *nbr,
                                 reads_info_t *reads_info)
{
        dpu_request_t *new_read =
                (dpu_request_t *) (reads->reads_area + reads->nb_reads * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes));
        new_read->offset = offset;
        new_read->count = count;
        new_read->num = num;

        memcpy(((uint8_t *)new_read) + sizeof(dpu_request_t), nbr, (size_t) reads_info->size_neighbour_in_bytes);
        reads->nb_reads++;
}

void add_seed_to_dpu_requests(dispatch_request_t *requests,
                              int num_read,
                              __attribute__((unused)) int nb_read_written,
                              index_seed_t *seed,
                              int8_t *nbr,
                              reads_info_t *reads_info)
{
        dispatch_request_t *this_request = requests + seed->num_dpu;
        dispatch_request_add(this_request,
                             (unsigned int) seed->offset,
                             (unsigned int) seed->nb_nbr,
                             (unsigned int) num_read,
                             nbr,
                             reads_info);
}

static void print_memory_layout(mram_info_t *mram_info, reads_info_t *reads_info)
{
        printf("\t                 addr       size\n");
        printf("\tmram_info        0x%.8x 0x%.8x\n", MRAM_INFO_ADDR, (unsigned int)sizeof(mram_info_t));
        printf("\tref inputs       0x%.8x 0x%.8x\n", (unsigned int)DPU_INPUTS_ADDR, mram_info->total_nbr_size);
        printf("\trequest_info     0x%.8x 0x%.8x\n",
               (unsigned int)DPU_REQUEST_INFO_ADDR(mram_info),
               (unsigned int)sizeof(request_info_t));
        printf("\trequest          0x%.8x 0x%.8x\n"
               "\t                 (one:  0x%.8x)\n",
               (unsigned int)DPU_REQUEST_ADDR(mram_info),
               (unsigned int)DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes) * MAX_DPU_REQUEST,
               (unsigned int)DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes));
        printf("\tempty space      0x%.8x 0x%.8x\n",
               (unsigned int)(DPU_REQUEST_ADDR(mram_info)
                              + DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes) * MAX_DPU_REQUEST),
               (unsigned int)(DPU_COMPUTE_TIME_ADDR
                              - (DPU_REQUEST_ADDR(mram_info)
                                 + DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes) * MAX_DPU_REQUEST)));
        printf("\tdpu time stats   0x%.8x 0x%.8x\n", (unsigned int)DPU_COMPUTE_TIME_ADDR,(unsigned int)DPU_COMPUTE_TIME_SIZE);
        printf("\ttasklet stats    0x%.8x 0x%.8x\n", (unsigned int)DPU_TASKLET_STATS_ADDR,(unsigned int)DPU_TASKLET_STATS_SIZE);
        printf("\tresult swap area 0x%.8x 0x%.8x\n", (unsigned int)DPU_SWAP_RESULT_ADDR, (unsigned int)DPU_SWAP_RESULT_SIZE);
        printf("\tresult area      0x%.8x 0x%.8x\n", (unsigned int)DPU_RESULT_ADDR, (unsigned int)DPU_RESULT_SIZE);

        assert((MRAM_INFO_ADDR + sizeof(mram_info_t)) <= DPU_INPUTS_ADDR);
        assert((DPU_INPUTS_ADDR + mram_info->total_nbr_size) <= DPU_REQUEST_INFO_ADDR(mram_info));
        assert((DPU_REQUEST_INFO_ADDR(mram_info) + sizeof(request_info_t)) <= DPU_REQUEST_ADDR(mram_info));
        assert((DPU_REQUEST_ADDR(mram_info) + MAX_DPU_REQUEST * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes))
               <= DPU_COMPUTE_TIME_ADDR);
        assert((DPU_COMPUTE_TIME_ADDR + DPU_COMPUTE_TIME_SIZE) <= DPU_TASKLET_STATS_ADDR);
        assert((DPU_TASKLET_STATS_ADDR + DPU_TASKLET_STATS_SIZE) <= DPU_SWAP_RESULT_ADDR);
        assert((DPU_SWAP_RESULT_ADDR + DPU_SWAP_RESULT_SIZE) <= DPU_RESULT_ADDR);
        assert((DPU_RESULT_SIZE + DPU_RESULT_ADDR) <= MRAM_SIZE);
}

void run_on_dpu(dispatch_request_t *dispatch,
                devices_t *devices,
                unsigned int dpu_offset,
                unsigned int nb_pass,
                sem_t *dispatch_free_sem,
                sem_t *acc_wait_sem,
                times_ctx_t *times_ctx,
                reads_info_t *reads_info)
{
        double t1, t2, t3, t4;
        unsigned int nb_ranks_per_run = devices->nb_ranks_per_run;
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;

        t1 = my_clock();

        PRINT_TIME_WRITE_READS(times_ctx, nb_pass);
        for (unsigned int each_rank = 0; each_rank < nb_ranks_per_run; each_rank++) {
                printf("() write inputs to rank #%d\n", each_rank);
                dpu_try_write_dispatch_into_mram(each_rank,
                                                 dpu_offset + each_rank * nb_dpus_per_rank,
                                                 devices,
                                                 dispatch,
                                                 reads_info);
        }
        sem_post(dispatch_free_sem);

        PRINT_TIME_WRITE_READS(times_ctx, nb_pass);
        t2 = my_clock();
        PRINT_TIME_COMPUTE(times_ctx, nb_pass);

        for (unsigned int each_rank = 0; each_rank < nb_ranks_per_run; each_rank++) {
                printf("() boot DPU rank #%d\n", each_rank);
                dpu_try_run(each_rank, devices);
        }

        for (unsigned int each_rank = 0; each_rank < nb_ranks_per_run; each_rank++) {
                while (!dpu_try_check_status(each_rank, devices));
                printf("DPU rank #%u completed\n", each_rank);
        }

        PRINT_TIME_COMPUTE(times_ctx, nb_pass);
        t3 = my_clock();
        PRINT_TIME_READ_RES(times_ctx, nb_pass);

        sem_wait(acc_wait_sem);
        /* Gather results and free DPUs */
        for (unsigned int each_rank = 0; each_rank < nb_ranks_per_run; each_rank++) {
                dpu_result_out_t *results[nb_dpus_per_rank];
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                        unsigned int this_dpu = each_rank * nb_dpus_per_rank + dpu_offset + each_dpu;
                        results[each_dpu] = get_mem_dpu_res(this_dpu);
                }
                dpu_try_get_results_and_log(each_rank, dpu_offset, devices, results);

                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                        unsigned int this_dpu = each_rank * nb_dpus_per_rank + dpu_offset + each_dpu;
                        if (dispatch[this_dpu].nb_reads == 0) {
                                results[each_dpu][0].num = -1;
                        }
                }
        }

        PRINT_TIME_READ_RES(times_ctx, nb_pass);
        t4 = my_clock();
        times_ctx->map_read = t4 - t1;
        times_ctx->tot_map_read += t4 - t1;
        times_ctx->write_reads = t2 - t1;
        times_ctx->tot_write_reads += t2 - t1;
        times_ctx->compute = t3 - t2;
        times_ctx->tot_compute += t3 - t2;
        times_ctx->read_result = t4 - t3;
        times_ctx->tot_read_result += t4 - t3;
}

void init_backend_dpu(devices_t **devices,
                             unsigned int nb_dpu_per_run,
                             const char *dpu_binary,
                             index_seed_t ***index_seed,
                             unsigned int nb_dpu,
                             __attribute__((unused)) genome_t *ref_genome,
                             __attribute__((unused)) reads_info_t *reads_info,
                             __attribute__((unused)) times_ctx_t *times_ctx,
                             __attribute__((unused)) backends_functions_t *backends_functions)
{
        malloc_dpu_res(nb_dpu);
        *index_seed = load_index_seeds();
        *devices = dpu_try_alloc_for(nb_dpu_per_run, dpu_binary);
}

void free_backend_dpu(devices_t *devices, __attribute__((unused)) unsigned int nb_dpu)
{
        dpu_try_free(devices);
        free_dpu_res();
}

void load_mram_dpu(unsigned int dpu_offset, devices_t *devices, reads_info_t *reads_info, times_ctx_t *times_ctx)
{
        double t1, t2;
        unsigned int nb_dpus_per_rank = devices->nb_dpus_per_rank;
        mram_info_t **mram = (mram_info_t **)malloc(nb_dpus_per_rank * sizeof(mram_info_t *));
        assert(mram != NULL);

        t1 = my_clock();

        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                mram[each_dpu] = (mram_info_t *)malloc(MRAM_SIZE);
                assert(mram[each_dpu] != NULL);
        }

        for (unsigned int each_rank = 0; each_rank < devices->nb_ranks_per_run; each_rank++) {
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                        unsigned int this_dpu = dpu_offset + each_dpu + each_rank * nb_dpus_per_rank;
                        printf("() load MRAM %d into DPU %d of rank %d\n", this_dpu, each_dpu, each_rank);
                        mram_load(mram[each_dpu], this_dpu);
                        print_memory_layout(mram[each_dpu], reads_info);
                }
                dpu_try_write_mram(each_rank, devices, mram);
        }
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_rank; each_dpu++) {
                free(mram[each_dpu]);
        }
        free(mram);

        t2 = my_clock();
        times_ctx->write_mram = t2 - t1;
        times_ctx->tot_write_mram += t2 - t1;
}
