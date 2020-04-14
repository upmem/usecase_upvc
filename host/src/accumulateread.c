#include "accumulateread.h"
#include "common.h"
#include "parse_args.h"
#include "upvc.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))

static FILE **result_file;
static dpu_result_out_t **result_list;
static acc_results_t *results_buffers[NB_DISPATCH_AND_ACC_BUFFER];
#define RESULTS_BUFFERS(pass_id) results_buffers[(pass_id) % NB_DISPATCH_AND_ACC_BUFFER]

SLIST_HEAD(min_heap_head, min_heap);
struct min_heap {
    dpu_result_out_t *list;
    SLIST_ENTRY(min_heap) entries;
};

static int cmpresult(const void *a, const void *b)
{
    dpu_result_out_t *A = (dpu_result_out_t *)a;
    dpu_result_out_t *B = (dpu_result_out_t *)b;
    if (A->num == B->num) {
        if (A->score > B->score) {
            return 1;
        } else {
            return -1;
        }
    }
    if (A->num > B->num) {
        return 1;
    } else {
        return -1;
    }
}

static void insert_in_order(struct min_heap_head *head, struct min_heap *elem)
{
    if (elem->list[0].num == -1)
        return;
    if (SLIST_EMPTY(head)) {
        SLIST_INSERT_HEAD(head, elem, entries);
        return;
    }
    struct min_heap *curr, *prev = NULL;
    SLIST_FOREACH (curr, head, entries) {
        if (cmpresult(&elem->list[0], &curr->list[0]) == -1) {
            if (prev == NULL)
                SLIST_INSERT_HEAD(head, elem, entries);
            else
                SLIST_INSERT_AFTER(prev, elem, entries);
            return;
        }
        prev = curr;
    }
    SLIST_INSERT_AFTER(prev, elem, entries);
}

static void init_min_heap_lists(
    struct min_heap_head *head, struct min_heap *min_inputs, unsigned int nb_input, dpu_result_out_t **inputs)
{
    SLIST_INIT(head);
    for (unsigned int each_input = 0; each_input < nb_input; each_input++) {
        min_inputs[each_input].list = inputs[each_input];
        insert_in_order(head, &min_inputs[each_input]);
    }
}

static void merge_result_list(
    dpu_result_out_t *output, dpu_result_out_t **inputs, unsigned int nb_input, unsigned int output_size)
{
    struct min_heap_head head;
    struct min_heap min_inputs[nb_input];
    init_min_heap_lists(&head, min_inputs, nb_input, inputs);
    for (unsigned int output_id = 0; output_id < output_size; output_id++) {
        assert(!SLIST_EMPTY(&head));
        struct min_heap *min_input = SLIST_FIRST(&head);
        SLIST_REMOVE_HEAD(&head, entries);

        output[output_id] = min_input->list[0];
        min_input->list = &min_input->list[1];

        insert_in_order(&head, min_input);
    }
}

acc_results_t accumulate_get_result(unsigned int pass_id)
{
    if (result_file[pass_id] == NULL) {
        static const dpu_result_out_t dummy_res = { .num = -1 };
        char result_filename[512];
        sprintf(result_filename, "result_%u.bin", pass_id);

        result_file[pass_id] = fopen(result_filename, "w+");
        fwrite(&dummy_res, sizeof(dummy_res), 1, result_file[pass_id]);
    }

    fseek(result_file[pass_id], 0, SEEK_END);
    size_t size = ftell(result_file[pass_id]);
    rewind(result_file[pass_id]);

    dpu_result_out_t *results = (dpu_result_out_t *)malloc(size);
    assert(results != NULL);
    size_t size_read = fread(results, size, 1, result_file[pass_id]);
    assert(size_read == 1);

    return (acc_results_t) { .nb_res = (size / sizeof(dpu_result_out_t)) - 1, .results = results };
}

void accumulate_read(unsigned int pass_id, unsigned int dpu_offset)
{
    unsigned int nb_dpus_per_run = backends_functions.get_nb_dpus_per_run();
    unsigned int nb_dpus = index_get_nb_dpu();
    acc_results_t *acc_res = RESULTS_BUFFERS(pass_id);

    unsigned int nb_dpus_used_current_run = MIN(nb_dpus - dpu_offset, nb_dpus_per_run);

    // sort the lists comming from each dpu independantly
#pragma omp parallel for
    for (unsigned int numdpu = 0; numdpu < nb_dpus_per_run; numdpu++) {
        if (numdpu + dpu_offset >= nb_dpus)
            continue;
        result_list[numdpu] = acc_res[numdpu].results;
        qsort(result_list[numdpu], acc_res[numdpu].nb_res, sizeof(dpu_result_out_t), cmpresult);
    }

    // compute the total number of resultat for all DPUs
    nb_result_t total_nb_res = 0;
    for (unsigned int numdpu = 0; numdpu < nb_dpus_per_run; numdpu++) {
        if (numdpu + dpu_offset >= nb_dpus)
            break;
        total_nb_res += acc_res[numdpu].nb_res;
        if (acc_res[numdpu].results[acc_res[numdpu].nb_res].num != -1) {
            ERROR_EXIT(-72, "%s:[P%u, M%u]: end mark is not there in DPU#%u\n", __func__, pass_id, dpu_offset, numdpu);
        }
    }

    // Get data from FILE *
    acc_results_t acc_res_from_file = accumulate_get_result(pass_id);
    unsigned int nb_read = acc_res_from_file.nb_res + total_nb_res;
    size_t size = sizeof(dpu_result_out_t) * (nb_read + 1);

    // Alloc the merged and sorted tab of result for this pass
    dpu_result_out_t *merged_result_tab = malloc(size);
    assert(merged_result_tab != NULL);

    // Merge and sort all the result for this pass
    result_list[nb_dpus_used_current_run] = acc_res_from_file.results;
    merge_result_list(merged_result_tab, result_list, nb_dpus_used_current_run + 1, nb_read);

    // update FILE *
    free(acc_res_from_file.results);
    merged_result_tab[nb_read].num = -1;
    rewind(result_file[pass_id]);
    size_t written_size = fwrite(merged_result_tab, size, 1, result_file[pass_id]);
    assert(written_size == 1);
    free(merged_result_tab);
}

void accumulate_free()
{
    unsigned int nb_dpu = backends_functions.get_nb_dpus_per_run();
    for (unsigned int each_pass = 0; each_pass < MAX_NB_PASS; each_pass++) {
        if (result_file[each_pass] != NULL)
            fclose(result_file[each_pass]);
    }
    free(result_file);

    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            free(results_buffers[each_pass][each_dpu].results);
        }
        free(results_buffers[each_pass]);
    }

    free(result_list);
}

void accumulate_init()
{
    unsigned int nb_dpu = backends_functions.get_nb_dpus_per_run();
    result_file = (FILE **)calloc(MAX_NB_PASS, sizeof(FILE *));
    assert(result_file != NULL);

    for (unsigned int each_pass = 0; each_pass < NB_DISPATCH_AND_ACC_BUFFER; each_pass++) {
        results_buffers[each_pass] = (acc_results_t *)malloc(sizeof(acc_results_t) * nb_dpu);
        assert(results_buffers[each_pass] != NULL);
        for (unsigned int each_dpu = 0; each_dpu < nb_dpu; each_dpu++) {
            results_buffers[each_pass][each_dpu].results = (dpu_result_out_t *)malloc(sizeof(dpu_result_out_t) * MAX_DPU_RESULTS);
            assert(results_buffers[each_pass][each_dpu].results != NULL);
        }
    }

    result_list = (dpu_result_out_t **)malloc(sizeof(dpu_result_out_t *) * (nb_dpu + 1));
    assert(result_list != NULL);
}

acc_results_t *accumulate_get_buffer(unsigned int dpu_id, unsigned int pass_id) { return &(RESULTS_BUFFERS(pass_id)[dpu_id]); }
