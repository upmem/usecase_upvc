#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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

index_seed_t **get_index_seed_dpu(unsigned int nb_dpu,
                                  __attribute__((unused)) genome_t *ref_genome,
                                  reads_info_t *reads_info,
                                  __attribute__((unused)) times_ctx_t *times_ctx,
                                  __attribute__((unused)) backends_functions_t *backends_funcions)
{
        malloc_dpu(reads_info, nb_dpu);
        return reload_mram_images_and_seeds(reads_info);
}

vmi_t *init_vmis_dpu(unsigned int nb_dpu)
{
        vmi_t *vmis = (vmi_t *) calloc(nb_dpu, sizeof(vmi_t));
        char *vmi_name = (char *) malloc(strlen("0000."));
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                (void) sprintf(vmi_name, "%04u", dpuno);
                vmi_create(vmi_name, vmis + dpuno);
        }
        free(vmi_name);
        return vmis;
}

static void dump_mdpu_images_into_mram_files(vmi_t *vmis, unsigned int *nbr_counts, unsigned int nb_dpu, reads_info_t *reads_info)
{
        mram_info_t *mram_image;
        printf("Creating MRAM images\n");
        mram_image = mram_create(reads_info);
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
                vmi_t *this_vmi = vmis + each_dpu;
                mram_reset(mram_image, reads_info);
                mram_copy_vmi(mram_image, this_vmi, nbr_counts[each_dpu], reads_info);
                mram_save(mram_image, each_dpu);
        }
        mram_free(mram_image);
}

void free_vmis_dpu(vmi_t *vmis, unsigned int nb_dpu, unsigned int *nb_neighbours, reads_info_t *reads_info)
{
        dump_mdpu_images_into_mram_files(vmis, nb_neighbours, nb_dpu, reads_info);
        for (unsigned int dpuno = 0; dpuno < nb_dpu; dpuno++) {
                vmi_delete(vmis + dpuno);
        }
        free(vmis);
}

void write_vmi_dpu(vmi_t *vmis, unsigned int dpuno, unsigned int k, int8_t *nbr, long coords, reads_info_t *reads_info)
{
        unsigned int size_neighbour_in_bytes = reads_info->size_neighbour_in_bytes;
        unsigned int out_len = ALIGN_DPU(sizeof(long) + size_neighbour_in_bytes);
        uint8_t *temp_buff = (uint8_t *) malloc(out_len);
        memset(temp_buff, 0, out_len);
        ((long *) temp_buff)[0] = coords;
        memcpy(temp_buff + sizeof(long), nbr, (size_t) size_neighbour_in_bytes);
        vmi_write(vmis + dpuno, k * out_len, temp_buff, out_len);
        free(temp_buff);
}

static void dispatch_request_add(dispatch_request_t *reads,
                                 unsigned int offset,
                                 unsigned int count,
                                 unsigned int num,
                                 int8_t *nbr,
                                 reads_info_t *reads_info)
{
        unsigned int new_read_list_size_in_bytes;
        dpu_request_t *new_read;
        unsigned int read_index = reads->nb_reads;
        reads->nb_reads++;

        if (reads->nb_reads >= MAX_NB_DPU_READ) {
                ERROR_EXIT(27, "*** buffer full - cannot register more reads on one DPU - aborting");
        }

        new_read_list_size_in_bytes = reads->nb_reads * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes);
        reads->reads_area = realloc(reads->reads_area, new_read_list_size_in_bytes);

        new_read = (dpu_request_t *) (reads->reads_area + read_index * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes));
        new_read->offset = offset;
        new_read->count = count;
        new_read->num = num;

        memcpy(((uint8_t *)new_read) + sizeof(dpu_request_t), nbr, (size_t) reads_info->size_neighbour_in_bytes);
}

void add_seed_to_dpu_requests(dispatch_request_t *requests,
                              int num_dpu,
                              int num_read,
                              __attribute__((unused)) int nb_read_written,
                              index_seed_t *seed,
                              int8_t *nbr,
                              reads_info_t *reads_info)
{
        dispatch_request_t *this_request = requests + num_dpu;
        dispatch_request_add(this_request,
                             (unsigned int) seed->offset,
                             (unsigned int) seed->nb_nbr,
                             (unsigned int) num_read,
                             nbr,
                             reads_info);
}

static void print_memory_layout(mram_info_t *mram_info, unsigned int nb_reads, reads_info_t *reads_info)
{
        printf("\t                  addr       size\n");
        printf("\tmram_info         0x%.8x 0x%.8x\n", MRAM_INFO_ADDR, (unsigned int)sizeof(mram_info_t));
        printf("\tinputs (ref nbrs) 0x%.8x 0x%.8x\n", (unsigned int)DPU_INPUTS_ADDR, mram_info->total_nbr_size);
        printf("\trequest_info      0x%.8x 0x%.8x\n",
               (unsigned int)DPU_REQUEST_INFO_ADDR(mram_info),
               (unsigned int)sizeof(request_info_t));
        printf("\trequest           0x%.8x 0x%.8x (0x%.8x)\n",
               (unsigned int)DPU_REQUEST_ADDR(mram_info),
               (unsigned int)DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes) * nb_reads,
               (unsigned int)DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes));
        printf("\ttasklet stats     0x%.8x 0x%.8x\n", (unsigned int)DPU_TASKLET_STATS_ADDR,(unsigned int)DPU_TASKLET_STATS_SIZE);
        printf("\tresult swap area  0x%.8x 0x%.8x\n", (unsigned int)DPU_SWAP_RESULT_ADDR, (unsigned int)DPU_SWAP_RESULT_SIZE);
        printf("\tresult area       0x%.8x 0x%.8x\n", (unsigned int)DPU_RESULT_ADDR, (unsigned int)DPU_RESULT_SIZE);

        assert((MRAM_INFO_ADDR + sizeof(mram_info_t)) <= DPU_INPUTS_ADDR);
        assert((DPU_INPUTS_ADDR + mram_info->total_nbr_size) <= DPU_REQUEST_INFO_ADDR(mram_info));
        assert((DPU_REQUEST_INFO_ADDR(mram_info) + sizeof(request_info_t)) <= DPU_REQUEST_ADDR(mram_info));
        assert((DPU_REQUEST_ADDR(mram_info) + nb_reads * DPU_REQUEST_SIZE(reads_info->size_neighbour_in_bytes))
               <= DPU_TASKLET_STATS_ADDR);
        assert((DPU_TASKLET_STATS_ADDR + DPU_TASKLET_STATS_SIZE) <= DPU_SWAP_RESULT_ADDR);
        assert((DPU_SWAP_RESULT_ADDR + DPU_SWAP_RESULT_SIZE) <= DPU_RESULT_ADDR);
        assert((DPU_RESULT_SIZE + DPU_RESULT_ADDR) <= MRAM_SIZE);
}

void run_on_dpu(dispatch_t dispatch,
                devices_t devices,
                unsigned int nb_dpu,
                times_ctx_t *times_ctx,
                reads_info_t *reads_info)
{
        double t1, t2;
        unsigned int nb_dpus_per_run = get_nb_dpus_per_run();

        t1 = my_clock();

        if (devices == NULL) {
                ERROR("Unable to alloc devices!");
        }

        for (unsigned int first_dpu = ((DEBUG_FIRST_RUN != -1) ? (DEBUG_FIRST_RUN * nb_dpus_per_run) : 0);
             first_dpu < nb_dpu;
             first_dpu += nb_dpus_per_run) {

                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                mram_info_t *mram;
                                printf("() write MRAM #%d %u reads\n", this_dpu, dispatch[this_dpu].nb_reads);
                                dpu_try_load_mram_number(this_dpu, each_dpu, devices, &mram, reads_info);
                                print_memory_layout(mram, dispatch[this_dpu].nb_reads, reads_info);
                                dpu_try_write_dispatch_into_mram(each_dpu, devices,
                                                                 dispatch[this_dpu].nb_reads,
                                                                 dispatch[this_dpu].reads_area,
                                                                 mram,
                                                                 reads_info);
                        }
                }

                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                printf("() boot DPU #%d\n", this_dpu);
                                dpu_try_run(each_dpu, devices);
                        }
                }

                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                while (!dpu_try_check_status(each_dpu, devices)) {
                                        sleep(1);
                                }
                                printf("DPU #%u completed\n", this_dpu);
                        }
                }

                /* Gather results and free DPUs */
                for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
                        unsigned int this_dpu = first_dpu + each_dpu;
                        if (dispatch[this_dpu].nb_reads != 0) {
                                dpu_try_log(each_dpu, devices);
                                dpu_result_out_t *results = dpu_try_get_results(each_dpu, devices);
                                int i;
                                for (i = 0; results[i].num != -1; i++) {
                                        if (i == MAX_DPU_RESULTS) {
                                                ERROR_EXIT(22, "BUG! no EOR marker found when parsing DPU results!\n");
                                        }

                                        long coords = ((long) results[i].seed_nr) |
                                                (((long) (results[i].seq_nr)) << 32);
                                        write_result(this_dpu, i, results[i].num, coords, results[i].score);
                                }
                                write_result(this_dpu, i, (unsigned int) -1, -1L, (unsigned int) -1);

                                free(results);
                        }
                }

        }

        t2 = my_clock();
        times_ctx->map_read = t2 - t1;
        times_ctx->tot_map_read += t2 - t1;
}
