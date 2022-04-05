/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "accumulateread.h"
#include "common.h"
#include "index.h"
#include "upvc.h"
#include "profiling.h"
#include "debug.h"

#include <assert.h>
#include <dpu_backend.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>
#include <unistd.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))

static FILE **result_file;
static acc_results_t *results_buffers[NB_DISPATCH_AND_ACC_BUFFER];
#define RESULTS_BUFFERS(pass_id) results_buffers[(pass_id) % NB_DISPATCH_AND_ACC_BUFFER]

#define BUCKET_SIZE (16)
#define NB_BUCKET (1 << BUCKET_SIZE)
#define BUCKET_MASK (NB_BUCKET - 1)
typedef struct bucket_elem {
    dpu_result_out_t *elem;
    struct bucket_elem *next;
} bucket_elem_t;

static bucket_elem_t **bucket_table_read, **bucket_table_write, **bucket_table_last;
static bucket_elem_t *bucket_elems;

#define ACCUMULATE_THREADS (8)
#define ACCUMULATE_THREADS_SLAVE (ACCUMULATE_THREADS - 1)
static pthread_barrier_t barrier;
static pthread_t thread_id[ACCUMULATE_THREADS_SLAVE];
static bool stop_threads = false;
static acc_results_t *acc_res;
static unsigned int nb_dpus_used_current_run;
static unsigned int *dpu_offset_res;
static unsigned int nb_pass;

static void merge_dpu_list_into_buckets(const unsigned int thread_id)
{
    for (unsigned int numdpu = thread_id; numdpu < nb_dpus_used_current_run; numdpu += ACCUMULATE_THREADS) {
        acc_results_t *curr_res = &acc_res[numdpu];
        for (unsigned int each_res = 0; each_res < curr_res->nb_res; each_res++) {
            bucket_elem_t *elem = &bucket_elems[each_res + dpu_offset_res[numdpu]];
            dpu_result_out_t *res = &(curr_res->results[each_res]);
            elem->elem = res;
            elem->next = __atomic_exchange_n(&bucket_table_read[res->key & BUCKET_MASK], elem, __ATOMIC_RELAXED);
        }
    }
}

static void *accumulate_read_thread_fct(void *arg)
{
    const unsigned int thread_id = (int)(uintptr_t)arg;

    pthread_barrier_wait(&barrier);
    while (!stop_threads) {
        merge_dpu_list_into_buckets(thread_id);
        pthread_barrier_wait(&barrier);
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

typedef void (*insert_bucket_fct_t)(bucket_elem_t *, unsigned int);

static void insert_bucket_first(bucket_elem_t *bucket, const unsigned int bucket_id)
{
    bucket->next = bucket_table_write[bucket_id];
    bucket_table_write[bucket_id] = bucket;
}

static void insert_bucket_last(bucket_elem_t *bucket, const unsigned int bucket_id)
{
    if (bucket_table_last[bucket_id] == NULL) {
        bucket_table_write[bucket_id] = bucket;
        bucket_table_last[bucket_id] = bucket;
    } else {
        bucket_table_last[bucket_id]->next = bucket;
        bucket_table_last[bucket_id] = bucket;
    }
    bucket->next = NULL;
}

static void do_bucket_sort(const unsigned int bucket_pass, insert_bucket_fct_t insert)
{
    for (int each_bucket = NB_BUCKET - 1; each_bucket >= 0; each_bucket--) {
        bucket_elem_t *bucket = bucket_table_read[each_bucket];
        while (bucket != NULL) {
            bucket_elem_t *next_bucket = bucket->next;
            insert(bucket, (bucket->elem->key >> (BUCKET_SIZE * bucket_pass)) & BUCKET_MASK);
            bucket = next_bucket;
        }
    }
}

static bucket_elem_t *get_next_read_in_bucket(bucket_elem_t *curr_elem, unsigned int *curr_bucket_id)
{
    unsigned int bucket_id = *curr_bucket_id;
    curr_elem = curr_elem->next;
    while (curr_elem == NULL) {
        curr_elem = bucket_table_read[bucket_id++];
        if (bucket_id > NB_BUCKET) {
            return NULL;
        }
    }
    *curr_bucket_id = bucket_id;
    return curr_elem;
}

static void copy_bucket_to_dest(dpu_result_out_t *dest, unsigned int nb_elem, bucket_elem_t *curr_elem, unsigned int *bucket_id)
{
    for (unsigned int each_res = 0; each_res < nb_elem; each_res++) {
        dest[each_res] = *(curr_elem->elem);
        curr_elem = get_next_read_in_bucket(curr_elem, bucket_id);
    }
}

static void merge_bucket_and_acc_list(
    dpu_result_out_t *dest, acc_results_t *acc_res, const unsigned nb_read, const unsigned int nb_read_in_bucket)
{
    unsigned int acc_res_idx = 0;
    unsigned int bucket_read_idx = 0;
    unsigned int bucket_id = 0;
    bucket_elem_t fake_first_elem = { .elem = NULL, .next = NULL };
    bucket_elem_t *bucket_elem = get_next_read_in_bucket(&fake_first_elem, &bucket_id);
    if (acc_res->nb_res == 0) {
        copy_bucket_to_dest(&dest[0], nb_read_in_bucket, bucket_elem, &bucket_id);
        return;
    }
    for (unsigned int each_read = 0; each_read < nb_read; each_read++) {
        dpu_result_out_t *src;
        if (acc_res->results[acc_res_idx].key > bucket_elem->elem->key) {
            src = bucket_elem->elem;
            if (++bucket_read_idx >= nb_read_in_bucket) {
                dest[each_read++] = *src;
                memcpy(
                    &dest[each_read], &acc_res->results[acc_res_idx], sizeof(dpu_result_out_t) * (acc_res->nb_res - acc_res_idx));
                return;
            }
            bucket_elem = get_next_read_in_bucket(bucket_elem, &bucket_id);
        } else {
            src = &acc_res->results[acc_res_idx];
            if (++acc_res_idx >= acc_res->nb_res) {
                dest[each_read++] = *src;
                copy_bucket_to_dest(&dest[each_read], nb_read_in_bucket - bucket_read_idx, bucket_elem, &bucket_id);
                return;
            }
        }
        dest[each_read] = *src;
    }
}

acc_results_t accumulate_get_result(unsigned int pass_id, bool free_results)
{
    if (result_file[pass_id] == NULL) {
        static const dpu_result_out_t dummy_res = { .num = -1 };
        char result_filename[512];
        sprintf(result_filename, "result_%u.bin", pass_id);

        result_file[pass_id] = fopen(result_filename, "w+");
        assert(result_file[pass_id] != NULL);
        assert(unlink(result_filename) == 0);
        fwrite(&dummy_res, sizeof(dummy_res), 1, result_file[pass_id]);
    }

    fseek(result_file[pass_id], 0, SEEK_END);
    size_t size = ftell(result_file[pass_id]);
    rewind(result_file[pass_id]);

    dpu_result_out_t *results = (dpu_result_out_t *)malloc(size);
    assert(results != NULL);
    size_t size_read = fread(results, size, 1, result_file[pass_id]);
    assert(size_read == 1);

    if (free_results) {
        fclose(result_file[pass_id]);
        result_file[pass_id] = NULL;
    }
    return (acc_results_t) { .nb_res = (size / sizeof(dpu_result_out_t)) - 1, .results = results };
}

void accumulate_read(unsigned int pass_id, unsigned int dpu_offset)
{
    printf("DPU_OFFSET: %u - PASS_ID: %u\n", dpu_offset, pass_id);
    PRINT_ALL_FUNCTION_STAT();
    nb_dpus_used_current_run = MIN(index_get_nb_dpu() - dpu_offset, nb_dpus_per_run);
    acc_res = RESULTS_BUFFERS(pass_id);

    // compute the total number of resultat for all DPUs
    nb_result_t total_nb_res = 0;
    for (unsigned int numdpu = 0; numdpu < nb_dpus_used_current_run; numdpu++) {
        dpu_offset_res[numdpu] = total_nb_res;
        total_nb_res += acc_res[numdpu].nb_res;
        uint32_t rank, ci, dpu;
        get_dpu_info(numdpu, &rank, &ci, &dpu);
        if (acc_res[numdpu].results[acc_res[numdpu].nb_res].num != -1) {
            ERROR_EXIT(ERR_ACC_END_MARK_MISSING, "%s:[P%u, M%u]: end mark is not there in DPU#%u (0x%x.%u.%u)\n", __func__,
                pass_id, dpu_offset, numdpu, rank, ci, dpu);
        }
    }

    if (total_nb_res == 0) {
        return;
    }

    bucket_elems = (bucket_elem_t *)malloc(sizeof(bucket_elem_t) * total_nb_res);
    assert(bucket_elems != NULL);

    // merge the list comming from the DPUs in parallel in the bucket_table_read
    memset(bucket_table_read, 0, sizeof(bucket_elem_t *) * NB_BUCKET);
    pthread_barrier_wait(&barrier);
    merge_dpu_list_into_buckets(ACCUMULATE_THREADS_SLAVE);
    pthread_barrier_wait(&barrier);

    /* We need to make 4 passes as the key is 64bits and the bucket is (1 << 16) (64/16=4).
     * But the first pass has been done in parallel when mergind the list from all DPUs.
     * 3 passes left to do.
     */
    const insert_bucket_fct_t insert_fcts[4] = { NULL, insert_bucket_last, insert_bucket_last, insert_bucket_first };
    for (unsigned int bucket_pass = 1; bucket_pass < 4; bucket_pass++) {
        memset(bucket_table_write, 0, sizeof(bucket_elem_t *) * NB_BUCKET);
        memset(bucket_table_last, 0, sizeof(bucket_elem_t *) * NB_BUCKET);

        do_bucket_sort(bucket_pass, insert_fcts[bucket_pass]);

        bucket_elem_t **tmp = bucket_table_write;
        bucket_table_write = bucket_table_read;
        bucket_table_read = tmp;
    }

    // Get data from FILE *
    acc_results_t acc_res_from_file = accumulate_get_result(pass_id, false);
    unsigned int nb_read = acc_res_from_file.nb_res + total_nb_res;
    size_t size = sizeof(dpu_result_out_t) * (nb_read + 1);

    // Alloc the merged and sorted tab of result for this pass
    LOG_INFO("allocating %lu for dpu results\n", size);
    dpu_result_out_t *merged_result_tab = malloc(size);
    assert(merged_result_tab != NULL);

    // Merge and sort all the result for this pass
    merge_bucket_and_acc_list(merged_result_tab, &acc_res_from_file, nb_read, total_nb_res);

    // update FILE *
    free(acc_res_from_file.results);
    merged_result_tab[nb_read].num = -1;
    rewind(result_file[pass_id]);
    size_t written_size = fwrite(merged_result_tab, size, 1, result_file[pass_id]);
    assert(written_size == 1);
    free(merged_result_tab);
    free(bucket_elems);
}

void accumulate_free()
{
    for (unsigned int each_pass = 0; each_pass < nb_pass; each_pass++) {
        if (result_file[each_pass] != NULL)
            fclose(result_file[each_pass]);
    }
    free(result_file);

    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
            free(results_buffers[each_pass][each_dpu].results);
        }
        free(results_buffers[each_pass]);
    }

    stop_threads = true;
    pthread_barrier_wait(&barrier);
    assert(pthread_barrier_destroy(&barrier) == 0);
    for (unsigned int each_thread = 0; each_thread < ACCUMULATE_THREADS_SLAVE; each_thread++) {
        assert(pthread_join(thread_id[each_thread], NULL) == 0);
    }

    free(dpu_offset_res);

    free(bucket_table_read);
    free(bucket_table_write);
    free(bucket_table_last);
}

void accumulate_init(unsigned int max_nb_pass)
{
    nb_pass = max_nb_pass;
    check_ulimit_n(nb_pass + 16);

    result_file = (FILE **)calloc(nb_pass, sizeof(FILE *));
    assert(result_file != NULL);

    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        results_buffers[each_pass] = (acc_results_t *)malloc(sizeof(acc_results_t) * nb_dpus_per_run);
        assert(results_buffers[each_pass] != NULL);
        for (unsigned int each_dpu = 0; each_dpu < nb_dpus_per_run; each_dpu++) {
            results_buffers[each_pass][each_dpu].results = (dpu_result_out_t *)malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
            assert(results_buffers[each_pass][each_dpu].results != NULL);
        }
    }

    dpu_offset_res = (unsigned int *)malloc(sizeof(unsigned int) * nb_dpus_per_run);
    assert(dpu_offset_res != NULL);

    bucket_table_read = (bucket_elem_t **)malloc(sizeof(bucket_elem_t *) * NB_BUCKET);
    assert(bucket_table_read != NULL);
    bucket_table_write = (bucket_elem_t **)malloc(sizeof(bucket_elem_t *) * NB_BUCKET);
    assert(bucket_table_write != NULL);
    bucket_table_last = (bucket_elem_t **)malloc(sizeof(bucket_elem_t *) * NB_BUCKET);
    assert(bucket_table_last != NULL);

    assert(pthread_barrier_init(&barrier, NULL, ACCUMULATE_THREADS) == 0);
    for (unsigned int each_thread = 0; each_thread < ACCUMULATE_THREADS_SLAVE; each_thread++) {
        assert(pthread_create(&thread_id[each_thread], NULL, accumulate_read_thread_fct, (void *)(uintptr_t)each_thread) == 0);
    }
}

acc_results_t *accumulate_get_buffer(unsigned int dpu_id, unsigned int pass_id) { return &(RESULTS_BUFFERS(pass_id)[dpu_id]); }
