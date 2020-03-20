/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <semaphore.h>

#include "common.h"
#include "dispatch.h"
#include "dpus_mgmt.h"
#include "index.h"
#include "parse_args.h"
#include "upvc.h"
#include "upvc_dpu.h"

static void dispatch_request_add(
    dispatch_request_t *reads, unsigned int offset, unsigned int count, unsigned int num, int8_t *nbr)
{
    dpu_request_t *new_read = (dpu_request_t *)(reads->reads_area + reads->nb_reads * sizeof(dpu_request_t));
    new_read->offset = offset;
    new_read->count = count;
    new_read->num = num;

    memcpy(&new_read->nbr[0], nbr, (size_t)SIZE_NEIGHBOUR_IN_BYTES);
    reads->nb_reads++;
}

void add_seed_to_dpu_requests(
    dispatch_request_t *requests, int num_read, __attribute__((unused)) int nb_read_written, index_seed_t *seed, int8_t *nbr)
{
    dispatch_request_t *this_request = requests + seed->num_dpu;
    dispatch_request_add(this_request, (unsigned int)seed->offset, (unsigned int)seed->nb_nbr, (unsigned int)num_read, nbr);
}

void run_on_dpu(unsigned int dpu_offset, unsigned int rank_id, __attribute__((unused)) int delta_neighbour,
    sem_t *dispatch_free_sem, sem_t *acc_wait_sem)
{
    double t1, t2, t3, t4, t5;
    PRINT_TIME_WRITE_READS(rank_id);
    t1 = my_clock();
    dpu_try_write_dispatch_into_mram(rank_id, dpu_offset);
    sem_post(dispatch_free_sem);

    t2 = my_clock();
    PRINT_TIME_WRITE_READS(rank_id);
    PRINT_TIME_COMPUTE(rank_id);

    dpu_try_run(rank_id);

    PRINT_TIME_COMPUTE(rank_id);
    t3 = my_clock();

    sem_wait(acc_wait_sem);

    t4 = my_clock();
    PRINT_TIME_READ_RES(rank_id);
    dpu_try_get_results_and_log(rank_id, dpu_offset);
    PRINT_TIME_READ_RES(rank_id);
    t5 = my_clock();

    print(rank_id, "W %.3lf", t2 - t1);
    print(rank_id, "C %.3lf", t3 - t2);
    print(rank_id, "R %.3lf", t5 - t4);
}

void init_backend_dpu()
{
    malloc_dpu_res(get_nb_dpu());
    dpu_try_alloc_for();
}

void free_backend_dpu(unsigned int nb_dpu)
{
    dpu_try_free();
    free_dpu_res(nb_dpu);
}

void load_mram_dpu(unsigned int dpu_offset, unsigned int rank_id, int delta_neighbour)
{
    dpu_try_write_mram(rank_id, dpu_offset, delta_neighbour);
}
