/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __UPVC_H__
#define __UPVC_H__

#define DEBUG_ROUND (-1)
#define DEBUG_PASS (-1)
#define DEBUG_DPU (-1)

#define VERSION           "VERSION 1.5"
#define SIZE_SEED         14
#define MAX_READS_BUFFER  1048576    /* Maximum number of read by round        */

#define COST_SUB         10
#define COST_GAPO        11
#define COST_GAPE        1
#define NB_DIAG          15

#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13

#define CODE_A    0    /* ('A'>>1)&3   41H  0100 0001 */
#define CODE_C    1    /* ('C'>>1)&3   43H  0100 0011 */
#define CODE_T    2    /* ('T'>>1)&3   54H  0101 0100 */
#define CODE_G    3    /* ('G'>>1)&3   47H  0100 0111 */

/**
 * @brief Structure with information on the size of the reads and its neighbour.
 */
typedef struct {
        int size_read;
        int size_neighbour_in_32bits_words;
        int size_neighbour_in_bytes;
        int delta_neighbour_in_bytes;
} reads_info_t;

#include <pthread.h>
#include <stdio.h>
/**
 * @brief Structure to get the time to do each part of the application.
 */
typedef struct {
        double get_genome;
        double index_genome;
        double get_reads;
        double dispatch_read;
        double map_read;
        double process_read;
        double acc_read;
        double write_mram;
        double write_reads;
        double read_result;
        double compute;
        double tot_get_reads;
        double tot_dispatch_read;
        double tot_map_read;
        double tot_process_read;
        double tot_acc_read;
        double tot_write_mram;
        double tot_write_reads;
        double tot_read_result;
        double tot_compute;
        double vcf;
        FILE *time_file;
        pthread_mutex_t time_file_mutex;
} times_ctx_t;

#include <stddef.h>
#include <sys/time.h>
static inline double my_clock(void)
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (1.0e-6 * t.tv_usec + t.tv_sec);
}

#define PRINT_TIME(times_ctx, str, pass)                                \
        do {                                                            \
                pthread_mutex_lock(&times_ctx->time_file_mutex);        \
                fprintf(times_ctx->time_file, str, my_clock(), pass);   \
                fflush(times_ctx->time_file);                           \
                pthread_mutex_unlock(&times_ctx->time_file_mutex);      \
        } while (0)

#define PRINT_TIME_GET_READS(times_ctx, pass)           \
        PRINT_TIME(times_ctx, "%lf, %lf\n", pass + 0.0)
#define PRINT_TIME_DISPATCH(times_ctx, pass)                    \
        PRINT_TIME(times_ctx, "%lf, , %lf\n", pass + 0.02)
#define PRINT_TIME_ACC_READ(times_ctx, pass)                    \
        PRINT_TIME(times_ctx, "%lf, , , %lf\n", pass + 0.04)
#define PRINT_TIME_PROCESS_READ(times_ctx, pass)                \
        PRINT_TIME(times_ctx, "%lf, , , , %lf\n", pass + 0.06)
#define PRINT_TIME_WRITE_MRAM(times_ctx, pass, rank)                    \
        PRINT_TIME(times_ctx, "%lf, , , , , %lf\n", pass + 0.08 + rank * 0.1)
#define PRINT_TIME_WRITE_READS(times_ctx, pass, rank)                   \
        PRINT_TIME(times_ctx, "%lf, , , , , , %lf\n", pass + 0.10 + rank * 0.1)
#define PRINT_TIME_COMPUTE(times_ctx, pass, rank)                       \
        PRINT_TIME(times_ctx, "%lf, , , , , , , %lf\n", pass + 0.12 + rank * 0.1)
#define PRINT_TIME_READ_RES(times_ctx, pass, rank)                      \
        PRINT_TIME(times_ctx, "%lf, , , , , , , , %lf\n", pass + 0.14 + rank * 0.1)
#define PRINT_TIME_MAP_READ(times_ctx, pass, rank)                      \
        PRINT_TIME(times_ctx, "%lf, , , , , , , , , %lf\n", pass + 0.16 + rank * 0.1)

#include <stdlib.h>

#define WARNING(fmt, ...) do { fprintf(stderr, "WARNING: " fmt "\n", ##__VA_ARGS__);} while(0)
#define ERROR(fmt, ...) do { fprintf(stderr, "ERROR: " fmt "\n", ##__VA_ARGS__);} while(0)
#define ERROR_EXIT(err_code, fmt, ...) do { ERROR(fmt, ##__VA_ARGS__); exit((err_code));} while(0)

#endif /* __UPVC_H__ */
