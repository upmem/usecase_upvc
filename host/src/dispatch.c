/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "backends_functions.h"
#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "parse_args.h"
#include "upvc.h"

static dispatch_request_t *requests_buffers[NB_DISPATCH_AND_ACC_BUFFER];
#define REQUESTS_BUFFERS(pass_id) requests_buffers[(pass_id) % NB_DISPATCH_AND_ACC_BUFFER]

void dispatch_init(unsigned int nb_dpu)
{

    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        requests_buffers[each_pass] = (dispatch_request_t *)calloc(nb_dpu, sizeof(dispatch_request_t));
        assert(requests_buffers[each_pass] != NULL);

        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            requests_buffers[each_pass][each_dpu].dpu_requests = malloc(sizeof(dpu_request_t) * MAX_DPU_REQUEST);
            assert(requests_buffers[each_pass][each_dpu].dpu_requests != NULL);
        }
    }
}

void dispatch_free(unsigned int nb_dpu)
{
    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            free(requests_buffers[each_pass][each_dpu].dpu_requests);
        }
        free(requests_buffers[each_pass]);
    }
}

static void add_seed_to_requests(dispatch_request_t *reads, int num_read, index_seed_t *seed, int8_t *read, unsigned int num_dpu)
{
    dpu_request_t *new_read = &reads->dpu_requests[reads->nb_reads];
    new_read->offset = seed->offset;
    new_read->count = seed->nb_nbr;
    new_read->num = num_read;

    index_copy_neighbour((int8_t *)new_read->nbr, read);
    if (++reads->nb_reads > MAX_DPU_REQUEST) {
        ERROR_EXIT(28, "\nDispatch: Buffer full (DPU %d)", num_dpu);
    }
}

static void write_mem_DPU(
    dispatch_request_t *requests, index_seed_t *seed, int8_t *read, int num_read, pthread_mutex_t *dispatch_mutex)
{
    while (seed != NULL) {
        unsigned int num_dpu = seed->num_dpu;

        pthread_mutex_lock(&dispatch_mutex[num_dpu]);
        add_seed_to_requests(&requests[num_dpu], num_read, seed, read, num_dpu);
        pthread_mutex_unlock(&dispatch_mutex[num_dpu]);

        seed = seed->next;
    }
}

#define DISPATCHING_THREAD (8)

typedef struct {
    int8_t *read_buffer;
    int nb_read;
    pthread_mutex_t *dispatch_mutex;
    dispatch_request_t *requests;
} dispatch_thread_common_arg_t;

typedef struct {
    int thread_id;
    dispatch_thread_common_arg_t *common_arg;
} dispatch_thread_arg_t;

static void *dispatch_read_thread_fct(void *arg)
{
    dispatch_thread_arg_t *args = (dispatch_thread_arg_t *)arg;
    int thread_id = args->thread_id;
    int8_t *read_buffer = args->common_arg->read_buffer;
    int nb_read = args->common_arg->nb_read;
    pthread_mutex_t *dispatch_mutex = args->common_arg->dispatch_mutex;
    dispatch_request_t *requests = args->common_arg->requests;

    int nb_read_per_thread = nb_read / DISPATCHING_THREAD;
    for (int num_read = nb_read_per_thread * thread_id;
         (num_read < nb_read) && (num_read < (nb_read_per_thread * (thread_id + 1))); num_read++) {
        int8_t *read = &read_buffer[num_read * SIZE_READ];
        index_seed_t *seed = index_get(read);
        write_mem_DPU(requests, seed, read, num_read, dispatch_mutex);
    }
    return NULL;
}

void dispatch_read(unsigned int pass_id)
{
    int ret;
    int nb_dpu = get_nb_dpu();
    pthread_mutex_t dispatch_mutex[nb_dpu];
    dispatch_request_t *requests = REQUESTS_BUFFERS(pass_id);
    dispatch_thread_common_arg_t common_arg = {
        .read_buffer = get_reads_buffer(pass_id),
        .nb_read = get_reads_in_buffer(pass_id),
        .dispatch_mutex = dispatch_mutex,
        .requests = requests,
    };
    pthread_t thread_id[DISPATCHING_THREAD];
    dispatch_thread_arg_t thread_arg[DISPATCHING_THREAD];

    for (int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        requests[numdpu].nb_reads = 0;
        pthread_mutex_init(&dispatch_mutex[numdpu], NULL);
    }

    for (unsigned int each_dispatch_thread = 0; each_dispatch_thread < DISPATCHING_THREAD; each_dispatch_thread++) {
        thread_arg[each_dispatch_thread].thread_id = each_dispatch_thread;
        thread_arg[each_dispatch_thread].common_arg = &common_arg;
        ret = pthread_create(&thread_id[each_dispatch_thread], NULL, dispatch_read_thread_fct, &thread_arg[each_dispatch_thread]);
        assert(ret == 0);
    }

    for (unsigned int each_dispatch_thread = 0; each_dispatch_thread < DISPATCHING_THREAD; each_dispatch_thread++) {
        ret = pthread_join(thread_id[each_dispatch_thread], NULL);
        assert(ret == 0);
    }

    for (int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        pthread_mutex_destroy(&dispatch_mutex[numdpu]);
    }
}

dispatch_request_t *dispatch_get(unsigned int dpu_id, unsigned int pass_id) { return &(REQUESTS_BUFFERS(pass_id)[dpu_id]); }
