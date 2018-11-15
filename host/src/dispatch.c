/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dispatch.h"
#include "index.h"
#include "upvc_dpu.h"
#include "code.h"
#include "upvc.h"

static void dispatch_request_reset(dispatch_request_t *reads)
{
        reads->nb_reads = 0;
        if (reads->reads_area != NULL) {
                free(reads->reads_area);
        }
}

static void dispatch_request_free(dispatch_request_t *reads)
{
        dispatch_request_reset(reads);
}

static dispatch_request_t *dispatch_create(unsigned int nb_dpu)
{
        dispatch_request_t *result = (dispatch_request_t *) calloc(nb_dpu, sizeof(dispatch_request_t));
        if (result == NULL) {
                ERROR_EXIT(26, "*** could not allocate memory to store %u dispatch requests - aborting!", nb_dpu);
        }

        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
                dispatch_request_reset(result + each_dpu);
        }
        return result;
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

static void add_seed_to_requests(dispatch_request_t *requests,
                                 int num_dpu,
                                 int num_read,
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

static void write_mem_DPU(index_seed_t *seed,
                          int *counted_read,
                          int8_t *read,
                          int num_read,
                          dispatch_request_t *requests,
                          reads_info_t *reads_info,
                          bool simulation_mode)
{
        int8_t buf[reads_info->size_neighbour_in_bytes];

        while (seed != NULL) {
                int nb_read_written = counted_read[seed->num_dpu];

                code_neighbour(&read[SIZE_SEED], buf, reads_info);

                if (simulation_mode) {
                        write_count(seed->num_dpu, nb_read_written, seed->nb_nbr);
                        write_offset(seed->num_dpu, nb_read_written, seed->offset);
                        write_num(seed->num_dpu, nb_read_written, num_read);
                        write_neighbour_read(seed->num_dpu, nb_read_written, buf, reads_info);
                } else {
                        add_seed_to_requests(requests, seed->num_dpu, num_read, seed, buf, reads_info);
                }

                counted_read[seed->num_dpu]++;

                if (counted_read[seed->num_dpu] >= MAX_NB_DPU_READ - 1) {
                        ERROR_EXIT(28, "\nBuffer full (DPU %d)", seed->num_dpu);
                }
                seed = seed->next;
        }
}

void dispatch_free(dispatch_t requests, unsigned int nb_dpu)
{
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
                dispatch_request_free(requests + each_dpu);
        }
        free(requests);
}

dispatch_t dispatch_read(index_seed_t **index_seed,
                         int8_t* read_buffer,
                         int nb_read,
                         int nb_dpu,
                         times_ctx_t *times_ctx,
                         reads_info_t *reads_info,
                         bool simulation_mode)
{
        double t1, t2;
        int counted_read[nb_dpu];
        int size_read = reads_info->size_read;
        dispatch_request_t *requests = NULL;

        t1 = my_clock();

        requests = dispatch_create(nb_dpu);

        memset(&counted_read, 0, nb_dpu * sizeof(unsigned int));

        for (int num_read = 0; num_read < nb_read; num_read++) {
                int8_t *read = &read_buffer[num_read * size_read];
                index_seed_t *seed = index_seed[code_seed(read)];
                write_mem_DPU(seed,
                              counted_read,
                              read,
                              num_read,
                              requests,
                              reads_info,
                              simulation_mode);
        }

        for (int numdpu = 0; numdpu < nb_dpu; numdpu++) {
                int nb_counted_read = counted_read[numdpu];
                write_num(numdpu, nb_counted_read, -1);
        }

        t2 = my_clock();
        times_ctx->dispatch_read = t2 - t1;
        times_ctx->tot_dispatch_read += t2 - t1;

        return requests;
}
