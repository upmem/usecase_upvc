/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "upvc.h"

static dispatch_request_t *requests_buffers[NB_DISPATCH_AND_ACC_BUFFER];
#define REQUESTS_BUFFERS(pass_id) requests_buffers[(pass_id) % NB_DISPATCH_AND_ACC_BUFFER]

static unsigned int dispatch_pass_id;
static int8_t *read_buffer;
static int nb_read;
static dispatch_request_t *requests;

void dispatch_init()
{
    unsigned int nb_dpu = index_get_nb_dpu();
    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        requests_buffers[each_pass] = (dispatch_request_t *)calloc(nb_dpu, sizeof(dispatch_request_t));
        assert(requests_buffers[each_pass] != NULL);

        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            requests_buffers[each_pass][each_dpu].dpu_requests = malloc(sizeof(dpu_request_t) * MAX_DPU_REQUEST);
            assert(requests_buffers[each_pass][each_dpu].dpu_requests != NULL);
        }
    }
}

void dispatch_free()
{
    unsigned int nb_dpu = index_get_nb_dpu();
    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            free(requests_buffers[each_pass][each_dpu].dpu_requests);
        }
        free(requests_buffers[each_pass]);
    }
}

static void write_mem_DPU(index_seed_t *seed, int8_t *read, int num_read)
{
    while (seed != NULL) {
        unsigned int num_dpu = seed->num_dpu;
        unsigned int nb_reads = __sync_fetch_and_add(&requests[num_dpu].nb_reads, 1);
        dpu_request_t *new_read = &requests[num_dpu].dpu_requests[nb_reads];
        new_read->offset = seed->offset;
        new_read->count = seed->nb_nbr;
        new_read->num = num_read;

        index_copy_neighbour((int8_t *)new_read->nbr, read);
        if (nb_reads > MAX_DPU_REQUEST) {
            ERROR_EXIT(28, "%s:[P%u]: Buffer full (DPU#%u)", __func__, dispatch_pass_id, num_dpu);
        }

        seed = seed->next;
    }
}

#define DISPATCHING_THREAD (8)

static void *dispatch_read_thread_fct(void *arg)
{
    int thread_id = *(int *)arg;

    int nb_read_per_thread = nb_read / DISPATCHING_THREAD;
    for (int num_read = nb_read_per_thread * thread_id;
         (num_read < nb_read) && (num_read < (nb_read_per_thread * (thread_id + 1))); num_read++) {
        int8_t *read = &read_buffer[num_read * SIZE_READ];
        index_seed_t *seed = index_get(read);
        write_mem_DPU(seed, read, num_read);
    }
    return NULL;
}

void dispatch_read(unsigned int pass_id)
{
    int ret;
    int nb_dpu = index_get_nb_dpu();
    pthread_t thread_id[DISPATCHING_THREAD];
    int thread_arg[DISPATCHING_THREAD];
    requests = REQUESTS_BUFFERS(pass_id);
    read_buffer = get_reads_buffer(pass_id);
    nb_read = get_reads_in_buffer(pass_id);
    dispatch_pass_id = pass_id;

    for (int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        requests[numdpu].nb_reads = 0;
    }

    for (unsigned int each_dispatch_thread = 0; each_dispatch_thread < DISPATCHING_THREAD; each_dispatch_thread++) {
        thread_arg[each_dispatch_thread] = each_dispatch_thread;
        ret = pthread_create(&thread_id[each_dispatch_thread], NULL, dispatch_read_thread_fct, &thread_arg[each_dispatch_thread]);
        assert(ret == 0);
    }

    for (unsigned int each_dispatch_thread = 0; each_dispatch_thread < DISPATCHING_THREAD; each_dispatch_thread++) {
        ret = pthread_join(thread_id[each_dispatch_thread], NULL);
        assert(ret == 0);
    }

}

dispatch_request_t *dispatch_get(unsigned int dpu_id, unsigned int pass_id) { return &(REQUESTS_BUFFERS(pass_id)[dpu_id]); }
