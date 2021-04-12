/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_H__
#define __UPVC_H__

#define VERSION "VERSION 1.8"
#define MAX_READS_BUFFER (256 * 1024) /* Maximum number of read by round        */
#define NB_READS_BUFFER (16) /* To be increase if enough legacy memory available */
#define NB_DISPATCH_AND_ACC_BUFFER (4) /* To be increase if enough legacy memory available */
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
#include <sys/resource.h>
#include <sys/time.h>

#include <stddef.h>
#include <time.h>
static inline double my_clock(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return (1.0e-9 * t.tv_nsec + t.tv_sec);
}

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

#define XSTR(s) STR(s)
#define STR(s) #s

#define NB_RANKS_MAX (64)
extern unsigned int nb_dpus_per_run;
#endif /* __UPVC_H__ */
