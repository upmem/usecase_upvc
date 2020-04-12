/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_H__
#define __UPVC_H__

#define VERSION "VERSION 1.7"
#define MAX_READS_BUFFER (1048576) /* Maximum number of read by round        */
#define MAX_NB_PASS (1600 * (1024 * 1024ULL) / MAX_READS_BUFFER) /* Whole genomee is less than 1600 pass of 1Mreq. */
#define NB_READS_BUFFER (16)
#define NB_DISPATCH_AND_ACC_BUFFER (4)
#define NB_ROUND (1)

#define COST_SUB 10
#define COST_GAPO 11
#define COST_GAPE 1
#define NB_DIAG 15

#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13
#define CODE_ERR 14

#define CODE_A 0 /* ('A'>>1)&3   41H  0100 0001 */
#define CODE_C 1 /* ('C'>>1)&3   43H  0100 0011 */
#define CODE_T 2 /* ('T'>>1)&3   54H  0101 0100 */
#define CODE_G 3 /* ('G'>>1)&3   47H  0100 0111 */

#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <stddef.h>
#include <time.h>
static inline double my_clock(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return (1.0e-9 * t.tv_nsec + t.tv_sec);
}

extern FILE *time_file;
extern pthread_mutex_t time_file_mutex;
#define PRINT_TIME(str, pass)                                                                                                    \
    do {                                                                                                                         \
        pthread_mutex_lock(&time_file_mutex);                                                                                    \
        fprintf(time_file, str, my_clock(), pass);                                                                               \
        fflush(time_file);                                                                                                       \
        pthread_mutex_unlock(&time_file_mutex);                                                                                  \
    } while (0)

#define PRINT_TIME_GET_READS() PRINT_TIME("%lf, %lf\n", 0.0)
#define PRINT_TIME_DISPATCH() PRINT_TIME("%lf, , %lf\n", 0.1)
#define PRINT_TIME_ACC_READ() PRINT_TIME("%lf, , , %lf\n", 0.2)
#define PRINT_TIME_PROCESS_READ() PRINT_TIME("%lf, , , , %lf\n", 0.3)
#define PRINT_TIME_WRITE_MRAM(rank) PRINT_TIME("%lf, , , , , %lf\n", 0.4 + rank * 0.1)
#define PRINT_TIME_WRITE_READS(rank) PRINT_TIME("%lf, , , , , , %lf\n", 0.4 + rank * 0.1)
#define PRINT_TIME_COMPUTE(rank) PRINT_TIME("%lf, , , , , , , %lf\n", 0.4 + rank * 0.1)
#define PRINT_TIME_READ_RES(rank) PRINT_TIME("%lf, , , , , , , , %lf\n", 0.4 + rank * 0.1)
#define PRINT_TIME_MAP_READ(rank) PRINT_TIME("%lf, , , , , , , , , %lf\n", 0.4 + rank * 0.1)

#include <stdlib.h>

#define WARNING(fmt, ...)                                                                                                        \
    do {                                                                                                                         \
        fprintf(stderr, "WARNING: " fmt "\n", ##__VA_ARGS__);                                                                    \
    } while (0)
#define ERROR(fmt, ...)                                                                                                          \
    do {                                                                                                                         \
        fprintf(stderr, "ERROR: " fmt "\n", ##__VA_ARGS__);                                                                      \
    } while (0)
#define ERROR_EXIT(err_code, fmt, ...)                                                                                           \
    do {                                                                                                                         \
        ERROR(fmt, ##__VA_ARGS__);                                                                                               \
        exit((err_code));                                                                                                        \
    } while (0)

#define XSTR(s) STR(s)
#define STR(s) #s

#define NB_RANKS_MAX (64)

#include "backends_functions.h"
extern backends_functions_t backends_functions;
#define TABULATION "          "
#define LINE       "----------"
#define SEPARATOR  '|'
#define ENDLINE    "|\n\0"
static inline void print(const uint32_t rank_id, const char *fmt, ...)
{
    uint32_t each_rank;
    va_list args;
    va_start(args, fmt);
    char str[NB_RANKS_MAX * (strlen(TABULATION) + strlen(XSTR(SEPARATOR))) + strlen(ENDLINE) + 1 /* for '\0' */];
    int str_i = 0;
    for (each_rank = 0; each_rank < rank_id; each_rank++) {
        str[str_i++] = SEPARATOR;
        memcpy(&str[str_i], TABULATION, strlen(TABULATION));
        str_i += strlen(TABULATION);
    }
    str[str_i++] = SEPARATOR;
    int fmt_size = vsprintf(&str[str_i], fmt, args);
    memcpy(&str[str_i + fmt_size], TABULATION, strlen(TABULATION) - fmt_size);
    str_i += strlen(TABULATION);
    for (each_rank++; each_rank < backends_functions.get_nb_ranks_per_run(); each_rank++) {
        str[str_i++] = SEPARATOR;
        memcpy(&str[str_i], TABULATION, strlen(TABULATION));
        str_i += strlen(TABULATION);
    }
    memcpy(&str[str_i], ENDLINE, strlen(ENDLINE) + 1 /* for '\0' */);
    fprintf(stdout, "%s", str);
}

static inline void print_line(const uint32_t rank_id) { print(rank_id, LINE); }

#endif /* __UPVC_H__ */
