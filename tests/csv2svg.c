#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/queue.h>

#define STRING_MAX_SIZE (1024)

typedef enum {
    type_get_reads = 0,
    type_dispatch = 1,
    type_accumulate_read = 2,
    type_process_read = 3,
    type_write_mram = 4,
    type_write_reads = 5,
    type_compute = 6,
    type_read_result = 7,
    type_map_read = 8,
    nb_type = 9,
} type_t;

typedef struct line {
    double x1;
    double y1;
    double x2;
    double y2;
    type_t type;
    SLIST_ENTRY(line) next;
} line_t;
SLIST_HEAD(line_head, line);

static void add_point2(struct line_head *line_list, struct line_head *first_point_head, line_t *first_point, type_t type,
    double time, double value, unsigned int *height, double *width_offset)
{
    first_point->x2 = time - *width_offset;
    first_point->y2 = value;
    first_point->type = type;
    SLIST_REMOVE(first_point_head, first_point, line, next);
    SLIST_INSERT_HEAD(line_list, first_point, next);

    if (first_point->y2 > (float)(*height)) {
        *height = (unsigned int)(first_point->y2 + 1);
    }
}

static void add_point1(struct line_head *first_point_head, double time, double value, double *width_offset)
{
    if (*width_offset == 0) {
        *width_offset = time;
    }
    line_t *new_point = (line_t *)malloc(sizeof(line_t));
    new_point->x1 = time - *width_offset;
    new_point->y1 = value;
    SLIST_INSERT_HEAD(first_point_head, new_point, next);
}

static void add_point(struct line_head *line_list, struct line_head *first_point_head, type_t type, double time, double value,
    unsigned int *height, double *width_offset)
{
    line_t *first_point;
    SLIST_FOREACH(first_point, &first_point_head[type], next)
    {
        if (value != first_point->y1) {
            continue;
        }
        add_point2(line_list, &first_point_head[type], first_point, type, time, value, height, width_offset);
        return;
    }
    add_point1(&first_point_head[type], time, value, width_offset);
}

static bool read_line(FILE *fp, char *str, struct line_head *line_list, struct line_head *first_point, unsigned int *width,
    unsigned int *height, double *width_offset)
{
    char *ret_str;
    double time, get_reads, dispatch, accumulate_read, process_read, write_mram, write_reads, compute, read_result, map_read;
    ret_str = fgets(str, STRING_MAX_SIZE, fp);

    if (ret_str == NULL) {
        return false;
    }

    if (sscanf(str, "%lf, %lf", &time, &get_reads) == 2) {
        add_point(line_list, first_point, type_get_reads, time, get_reads, height, width_offset);
    } else if (sscanf(str, "%lf, , %lf", &time, &dispatch) == 2) {
        add_point(line_list, first_point, type_dispatch, time, dispatch, height, width_offset);
    } else if (sscanf(str, "%lf, , , %lf", &time, &accumulate_read) == 2) {
        add_point(line_list, first_point, type_accumulate_read, time, accumulate_read, height, width_offset);
    } else if (sscanf(str, "%lf, , , , %lf", &time, &process_read) == 2) {
        add_point(line_list, first_point, type_process_read, time, process_read, height, width_offset);
    } else if (sscanf(str, "%lf, , , , , %lf", &time, &write_mram) == 2) {
        add_point(line_list, first_point, type_write_mram, time, write_mram, height, width_offset);
    } else if (sscanf(str, "%lf, , , , , , %lf", &time, &write_reads) == 2) {
        add_point(line_list, first_point, type_write_reads, time, write_reads, height, width_offset);
    } else if (sscanf(str, "%lf, , , , , , , %lf", &time, &compute) == 2) {
        add_point(line_list, first_point, type_compute, time, compute, height, width_offset);
    } else if (sscanf(str, "%lf, , , , , , , , %lf", &time, &read_result) == 2) {
        add_point(line_list, first_point, type_read_result, time, read_result, height, width_offset);
    } else if (sscanf(str, "%lf, , , , , , , , , %lf", &time, &map_read) == 2) {
        add_point(line_list, first_point, type_map_read, time, map_read, height, width_offset);
    } else {
        return false;
    }

    *width = (unsigned int)(time - *width_offset + 1);

    return true;
}

#define ZOOM (100)

static void generate_svg(FILE *fp, char *str, struct line_head *line_list, unsigned int width, unsigned int height)
{
    char *type2string[nb_type] = {
        [type_get_reads] = "red",
        [type_dispatch] = "blue",
        [type_accumulate_read] = "green",
        [type_process_read] = "brown",
        [type_write_mram] = "yellow",
        [type_write_reads] = "orange",
        [type_compute] = "purple",
        [type_read_result] = "grey",
        [type_map_read] = "pink",
    };
    double acc[nb_type] = { 0.0 };
    unsigned int nb_it[nb_type] = { 0 };
    line_t *current_line;
    sprintf(str, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    fwrite(str, sizeof(char), strlen(str), fp);
    sprintf(str,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" width=\"%u\" "
        "height=\"%u\">\n",
        width * ZOOM, height * ZOOM);
    fwrite(str, sizeof(char), strlen(str), fp);

    SLIST_FOREACH(current_line, line_list, next)
    {
        acc[current_line->type] += (current_line->x2 - current_line->x1);
        nb_it[current_line->type]++;
        sprintf(str, "\t<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" stroke=\"%s\" stroke-width=\"%lf\"/>\n",
            current_line->x1 * ZOOM, current_line->y1 * ZOOM, current_line->x2 * ZOOM, current_line->y2 * ZOOM,
            type2string[current_line->type], 0.01 * ZOOM);
        fwrite(str, sizeof(char), strlen(str), fp);
    }

    printf("get reads: %lfs\n", acc[type_get_reads] / nb_it[type_get_reads]);
    printf("dispatch: %lfs\n", acc[type_dispatch] / nb_it[type_dispatch]);
    printf("accumulate_read: %lfs\n", acc[type_accumulate_read] / nb_it[type_accumulate_read]);
    printf("process_read: %lfs\n", acc[type_process_read] / nb_it[type_process_read]);
    printf("map read: %lfs\n", acc[type_map_read] / nb_it[type_map_read]);

    sprintf(str, "</svg>");
    fwrite(str, sizeof(char), strlen(str), fp);
}

int main(__attribute((unused)) int argc, char **argv)
{
    char str[STRING_MAX_SIZE];
    char *ret_str;
    struct line_head first_point[nb_type];
    struct line_head line_list = SLIST_HEAD_INITIALIZER(line_list);
    unsigned int width = 0, height = 0;
    double width_offset = 0.0;

    for (unsigned int each_type = 0; each_type < nb_type; each_type++) {
        SLIST_INIT(&first_point[each_type]);
    }

    {
        FILE *fp_in = fopen(argv[1], "r");

        /* Read first line */
        ret_str = fgets(str, STRING_MAX_SIZE, fp_in);
        assert(ret_str != NULL);

        while (read_line(fp_in, str, &line_list, first_point, &width, &height, &width_offset))
            ;
        fclose(fp_in);
    }
    printf("read complete!\n");

    {
        FILE *fp_out;
        sprintf(str, "%s.svg", argv[1]);
        fp_out = fopen(str, "w");

        generate_svg(fp_out, str, &line_list, width, height);

        fclose(fp_out);
    }
    printf("generate svg complete\n");

    return 0;
}
