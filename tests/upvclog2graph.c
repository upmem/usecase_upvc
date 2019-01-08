#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/queue.h>

#define STRING_MAX_SIZE (1024)

typedef struct dpu {
        unsigned long long cycle;
        unsigned int req;
        unsigned int nodp;
        unsigned int odpd;
        unsigned int results;
        unsigned int data_in;
        unsigned int result_out;
        unsigned int load;
        unsigned int store;
        SLIST_ENTRY(dpu) next;
} dpu_t;
SLIST_HEAD(dpu_head, dpu);

static bool read_new_line(FILE *fp, char *str, bool *end_of_file)
{
        char *ret_str = fgets(str, STRING_MAX_SIZE, fp);
        if (ret_str == NULL) {
                *end_of_file = true;
                return false;
        }
        return true;
}

#define LOG_PREFIX "LOG DPU="
#define GET_DPU_ID_FORMAT LOG_PREFIX "%u"
#define TID_FORMAT " TID=%u"
static bool read_dpu(FILE *fp, char *str, struct dpu_head *rank, bool *end_of_file)
{
        unsigned int curr_dpu_id;
        unsigned int dpu_id, tid, req, nodp, odpd, results, data_in, result_out, load, store;
        unsigned long long cycle;
        float time;
        dpu_t *new_dpu;

        if (*end_of_file) {
                return false;
        }
        if (sscanf(str, GET_DPU_ID_FORMAT, &curr_dpu_id) != 1) {
                return false;
        }

        new_dpu = (dpu_t *)calloc(1, sizeof(dpu_t));
        SLIST_INSERT_HEAD(rank, new_dpu, next);
        do {
                if (sscanf(str, GET_DPU_ID_FORMAT" TIME=%llu SEC=%f", &dpu_id, &cycle, &time) == 3) {
                        new_dpu->cycle = cycle;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " REQ=%u", &dpu_id, &tid, &req) == 3) {
                        new_dpu->req += req;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " NODP=%u", &dpu_id, &tid, &nodp) == 3) {
                        new_dpu->nodp += nodp;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " ODPD=%u", &dpu_id, &tid, &odpd) == 3) {
                        new_dpu->odpd += odpd;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " RESULTS=%u", &dpu_id, &tid, &results) == 3) {
                        new_dpu->results += results;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " DATA_IN=%u", &dpu_id, &tid, &data_in) == 3) {
                        new_dpu->data_in += data_in;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " RESULT_OUT=%u", &dpu_id, &tid, &result_out) == 3) {
                        new_dpu->result_out += result_out;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " LOAD=%u", &dpu_id, &tid, &load) == 3) {
                        new_dpu->load += load;
                } else if (sscanf(str, GET_DPU_ID_FORMAT TID_FORMAT " STORE=%u", &dpu_id, &tid, &store) == 3) {
                        new_dpu->store += store;
                } else {
                        return false;
                }
                if (!read_new_line(fp, str, end_of_file)) {
                        return false;
                }
                if (sscanf(str, GET_DPU_ID_FORMAT, &dpu_id) != 1) {
                        return false;
                }
        } while (dpu_id == curr_dpu_id);
        return true;
}

static bool read_unused(FILE *fp, char *str, bool *end_of_file)
{
        if (*end_of_file) {
                return false;
        }
        if (strncmp(str, LOG_PREFIX, strlen(LOG_PREFIX)) == 0) {
                return false;
        }
        return read_new_line(fp, str, end_of_file);
}

static void free_rank(struct dpu_head *rank)
{
        while (!SLIST_EMPTY(rank)) {
                dpu_t *dpu = SLIST_FIRST(rank);
                SLIST_REMOVE_HEAD(rank, next);
                free(dpu);
        }
}

#define EXPORT_MAX

static void export_rank(FILE *fp, char *str, struct dpu_head *rank)
{
        dpu_t *dpu;
#if defined EXPORT_MAX || EXPORT_MOY
        unsigned long long cycle = 0ULL;
        unsigned long long req = 0ULL, nodp = 0ULL, odpd = 0ULL, results = 0ULL, data_in = 0ULL, result_out = 0ULL, load = 0ULL, store = 0ULL;
#ifdef EXPORT_MOY
        unsigned int nb_dpu = 0;
#endif
#endif
        SLIST_FOREACH(dpu, rank, next) {
#ifdef EXPORT_ALL
                sprintf(str, "%llu, %u, %u, %u, %u, %u, %u, %u, %u\n",
                        dpu->cycle,
                        dpu->req,
                        dpu->nodp,
                        dpu->odpd,
                        dpu->results,
                        dpu->data_in,
                        dpu->result_out,
                        dpu->load,
                        dpu->store);
                fwrite(str, sizeof(char), strlen(str), fp);
#else
#if defined EXPORT_MAX || EXPORT_MOY
#ifdef EXPORT_MAX
#define FCT(a, b) a = b > a ? b : a
#elif EXPORT_MOY
#define FCT(a, b) a += b
                nb_dpu++;
#endif
                FCT(cycle, dpu->cycle);
                FCT(req, dpu->req);
                FCT(nodp, dpu->nodp);
                FCT(odpd, dpu->odpd);
                FCT(results, dpu->results);
                FCT(data_in, dpu->data_in);
                FCT(result_out, dpu->result_out);
                FCT(load, dpu->load);
                FCT(store, dpu->store);
#endif
#endif
        }
#ifdef EXPORT_MAX
        sprintf(str, "%llu, %llu, %llu, %llu, %llu, %llu, %llu, %llu, %llu\n",
                cycle,
                req,
                nodp,
                odpd,
                results,
                data_in,
                result_out,
                load,
                store);
        fwrite(str, sizeof(char), strlen(str), fp);
#elif EXPORT_MOY
        sprintf(str, "%llu, %llu, %llu, %llu, %llu, %llu, %llu, %llu, %llu\n",
                cycle / nb_dpu,
                req / nb_dpu,
                nodp / nb_dpu,
                odpd / nb_dpu,
                results / nb_dpu,
                data_in / nb_dpu,
                result_out / nb_dpu,
                load / nb_dpu,
                store / nb_dpu);
        fwrite(str, sizeof(char), strlen(str), fp);
#endif
}

int main(__attribute((unused)) int argc, char **argv)
{
        char str[STRING_MAX_SIZE];
        FILE *fp_in = fopen(argv[1], "r");
        FILE *fp_out;
        bool end_of_file = false;

        if (!read_new_line(fp_in, str, &end_of_file)) {
                return -1;
        }

        sprintf(str, "%s.csv", argv[1]);
        fp_out = fopen(str, "w");

        sprintf(str, "cycle, req, nodp, odpd, results, data_in, result_out, load, store\n");
        fwrite(str, sizeof(char), strlen(str), fp_out);

        while (!end_of_file) {
                struct dpu_head new_rank = SLIST_HEAD_INITIALIZER(new_rank);

                while (read_unused(fp_in, str, &end_of_file));
                while (read_dpu(fp_in, str, &new_rank, &end_of_file));

                export_rank(fp_out, str, &new_rank);
                free_rank(&new_rank);
        }
        fclose(fp_in);
        fclose(fp_out);

        printf("complete!\n");

        return 0;
}
