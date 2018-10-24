#include <stdio.h>
#include <stdlib.h>

#include "index.h"
#include "upvc_dpu.h"
#include "code.h"
#include "upvc.h"

static void write_mem_DPU(index_seed_t *seed, int *counted_read, int8_t *read, int num_read, reads_info_t *reads_info)
{
        int8_t buf[reads_info->size_neighbour_in_bytes];

        while (seed != NULL) {
                int nb_read_written = counted_read[seed->num_dpu];

                // TODO: to be replace by tranfer with DPUs memory
                write_count(seed->num_dpu, nb_read_written, seed->nb_nbr);
                write_offset(seed->num_dpu, nb_read_written, seed->offset);
                write_num(seed->num_dpu, nb_read_written, num_read);
                code_neighbour(&read[SIZE_SEED], buf, reads_info);
                write_neighbour_read(seed->num_dpu, nb_read_written, buf, reads_info);

                counted_read[seed->num_dpu]++;

                if (counted_read[seed->num_dpu] >= MAX_NB_DPU_READ - 1) {
                        printf("\nBuffer full (DPU %d)\n", seed->num_dpu);
                        exit(255);
                }
                seed = seed->next;
        }
}

void dispatch_read(index_seed_t **index_seed,
                   int8_t* read_buffer,
                   int nb_read,
                   times_ctx_t *times_ctx,
                   reads_info_t *reads_info)
{
        double t1, t2;
        int counted_read[NB_DPU] = {0};
        int size_read = reads_info->size_read;

        t1 = my_clock();

        for (int num_read = 0; num_read < nb_read; num_read++) {
                int8_t *read = &read_buffer[num_read * size_read];
                index_seed_t *seed = index_seed[code_seed(read)];
                write_mem_DPU(seed, counted_read, read, num_read, reads_info);
        }

        for (int numdpu = 0; numdpu < NB_DPU; numdpu++) {
                int nb_counted_read = counted_read[numdpu];
                write_num(numdpu, nb_counted_read, -1);
        }

        t2 = my_clock();
        times_ctx->dispatch_read = t2 - t1;
        times_ctx->tot_dispatch_read += t2 - t1;
}
