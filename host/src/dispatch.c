/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "backends_functions.h"
#include "code.h"
#include "dispatch.h"
#include "getread.h"
#include "index.h"
#include "parse_args.h"
#include "upvc.h"

#define DISPATCHING_THREAD (8)

static dispatch_request_t *requests;

void dispatch_init(unsigned int nb_dpu)
{
    requests = (dispatch_request_t *)calloc(nb_dpu, sizeof(dispatch_request_t));
    assert(requests != NULL);

    for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
        requests[each_dpu].reads_area = malloc(sizeof(dpu_request_t) * MAX_DPU_REQUEST);
    }
}

void dispatch_free(unsigned int nb_dpu)
{
    for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
        free(requests[each_dpu].reads_area);
    }
    free(requests);
}

static void write_mem_DPU(index_seed_t *seed, int *counted_read, int8_t *read, int num_read, pthread_mutex_t *dispatch_mutex)
{
    int8_t buf[SIZE_NEIGHBOUR_IN_BYTES];

    while (seed != NULL) {
        pthread_mutex_lock(&dispatch_mutex[seed->num_dpu]);
        {
            int nb_read_written = counted_read[seed->num_dpu];
            if ((nb_read_written + 1) > MAX_DPU_REQUEST) {
                ERROR_EXIT(28, "\nBuffer full (DPU %d)", seed->num_dpu);
            }

            code_neighbour(&read[SIZE_SEED], buf);
            backends_functions.add_seed_to_requests(requests, num_read, nb_read_written, seed, buf);
            counted_read[seed->num_dpu]++;
        }
        pthread_mutex_unlock(&dispatch_mutex[seed->num_dpu]);

        seed = seed->next;
    }
}

typedef struct {
    index_seed_t **index_seed;
    int8_t *read_buffer;
    int nb_read;
    int *counted_read;
    pthread_mutex_t *dispatch_mutex;
} dispatch_thread_common_arg_t;

typedef struct {
    int thread_id;
    dispatch_thread_common_arg_t *common_arg;
} dispatch_thread_arg_t;

static void *dispatch_read_thread_fct(void *arg)
{
    dispatch_thread_arg_t *args = (dispatch_thread_arg_t *)arg;
    int thread_id = args->thread_id;
    index_seed_t **index_seed = args->common_arg->index_seed;
    int8_t *read_buffer = args->common_arg->read_buffer;
    int nb_read = args->common_arg->nb_read;
    int *counted_read = args->common_arg->counted_read;
    pthread_mutex_t *dispatch_mutex = args->common_arg->dispatch_mutex;

    int nb_read_per_thread = nb_read / DISPATCHING_THREAD;
    for (int num_read = nb_read_per_thread * thread_id;
         (num_read < nb_read) && (num_read < (nb_read_per_thread * (thread_id + 1))); num_read++) {
        int8_t *read = &read_buffer[num_read * SIZE_READ];
        index_seed_t *seed = index_seed[code_seed(read)];
        write_mem_DPU(seed, counted_read, read, num_read, dispatch_mutex);
    }
    return NULL;
}

void dispatch_read(unsigned int pass_id)
{
    int ret;
    int nb_dpu = get_nb_dpu();
    int counted_read[nb_dpu];
    pthread_mutex_t dispatch_mutex[nb_dpu];
    dispatch_thread_common_arg_t common_arg = {
        .index_seed = index_get(),
        .read_buffer = get_reads_buffer(pass_id),
        .nb_read = get_reads_in_buffer(pass_id),
        .counted_read = counted_read,
        .dispatch_mutex = dispatch_mutex,
    };
    pthread_t thread_id[DISPATCHING_THREAD];
    dispatch_thread_arg_t thread_arg[DISPATCHING_THREAD];

    for (int numdpu = 0; numdpu < nb_dpu; numdpu++) {
        counted_read[numdpu] = 0;
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

dispatch_request_t *dispatch_get(unsigned int dpu_id) { return &requests[dpu_id]; }
