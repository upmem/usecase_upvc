/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_H__
#define __UPVC_H__

#define VERSION "VERSION 1.8"
#define MAX_READS_BUFFER (1048576) /* Maximum number of read by round        */
#define NB_READS_BUFFER (16)
#define NB_DISPATCH_AND_ACC_BUFFER (4)
#define NB_ROUND (1)

#define COST_SUB 10
#define COST_GAPO 11
#define COST_GAPE 1
#define NB_DIAG 15

#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

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

enum error_code {
    ERR_USAGE = -1,
    ERR_ACC_END_MARK_MISSING = -2,
    ERR_DISPATCH_BUFFER_FULL = -3,
    ERR_INDEX_FOLDER_ALREAD_EXIST = -4,
    ERR_INDEX_FOLDER_MKDIR_FAILED = -5,
    ERR_PROCESSREAD_DPD_FAILED = -6,
    ERR_SIMU_MAX_RESULTS_REACHED = -7,
    ERR_NO_GOAL_DEFINED = -8,
    ERR_CURRENT_FOLDER_PERMISSIONS = -8,
    ERR_FOPEN_FAILED = -9,
};

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

#define CHECK_FILE(f, name)                                                                                                      \
    do {                                                                                                                         \
        if (f == NULL)                                                                                                           \
            ERROR_EXIT(ERR_FOPEN_FAILED, "Could not open file '%s' (%s)", name, strerror(errno));                                \
    } while (0)

static inline void check_ulimit_n(unsigned int expected_limit)
{
    struct rlimit nofile_limit;
    assert(getrlimit(RLIMIT_NOFILE, &nofile_limit) == 0);
    if ((unsigned int)nofile_limit.rlim_cur < expected_limit) {
        WARNING("Number of file descriptor that can be opened by this process looks too small (current: %u - expected: %u), use "
                "'ulimit -n' to set to appropriate value",
            (unsigned int)nofile_limit.rlim_cur, expected_limit);
        printf("Press any key to continue\n");
        getchar();
    }
}

#define XSTR(s) STR(s)
#define STR(s) #s

#define NB_RANKS_MAX (64)

extern unsigned int nb_dpus_per_run;
extern unsigned int nb_ranks_per_run;
#define TABULATION "         "
#define LINE       "---------"
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
    for (each_rank++; each_rank < nb_ranks_per_run; each_rank++) {
        str[str_i++] = SEPARATOR;
        memcpy(&str[str_i], TABULATION, strlen(TABULATION));
        str_i += strlen(TABULATION);
    }
    memcpy(&str[str_i], ENDLINE, strlen(ENDLINE) + 1 /* for '\0' */);
    fprintf(stdout, "%s", str);
}

static inline void print_line(const uint32_t rank_id) { print(rank_id, LINE); }

#endif /* __UPVC_H__ */
