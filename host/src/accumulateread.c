#include "common.h"
#include "parse_args.h"
#include "upvc.h"
#include "upvc_dpu.h"

#include <sys/queue.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))

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

int accumulate_read(
    dpu_result_out_t **result_tab, unsigned int *result_tab_nb_read, unsigned int dpu_offset, times_ctx_t *times_ctx)
{
    double t1, t2;

    t1 = my_clock();
    unsigned int nb_dpus_per_run = get_nb_dpus_per_run();
    unsigned int nb_dpus = get_nb_dpu();
    nb_result_t *nb_res = get_mem_dpu_nb_res(dpu_offset);
    dpu_result_out_t **result_list = (dpu_result_out_t **)malloc(sizeof(dpu_result_out_t *) * (nb_dpus_per_run + 1));
    assert(result_list != NULL);
    unsigned int nb_dpus_used_current_run = MIN(get_nb_dpu() - dpu_offset, nb_dpus_per_run);
#pragma omp parallel for
    for (unsigned int numdpu = dpu_offset; numdpu < dpu_offset + nb_dpus_per_run; numdpu++) {
        if (numdpu >= nb_dpus)
            continue;
        dpu_result_out_t *dpu_res = get_mem_dpu_res(numdpu);
        qsort(dpu_res, nb_res[numdpu - dpu_offset], sizeof(dpu_result_out_t), cmpresult);
        result_list[numdpu - dpu_offset] = dpu_res;
    }
    nb_result_t total_nb_res = 0;
    for (unsigned int numdpu = dpu_offset; numdpu < dpu_offset + get_nb_dpus_per_run() && numdpu < get_nb_dpu(); numdpu++) {
        total_nb_res += nb_res[numdpu - dpu_offset];
    }
    *result_tab_nb_read += total_nb_res;
    dpu_result_out_t *merged_result_tab = malloc((*result_tab_nb_read + 1) * sizeof(dpu_result_out_t));
    assert(merged_result_tab != NULL);
    result_list[nb_dpus_used_current_run] = *result_tab;

    merge_result_list(merged_result_tab, result_list, nb_dpus_used_current_run + 1, *result_tab_nb_read);

    free(*result_tab);
    *result_tab = merged_result_tab;
    merged_result_tab[*result_tab_nb_read].num = -1;

    t2 = my_clock();
    times_ctx->acc_read = t2 - t1;
    times_ctx->tot_acc_read += t2 - t1;

    return (int)total_nb_res;
}
