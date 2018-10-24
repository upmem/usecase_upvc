#ifndef __UPVC_H__
#define __UPVC_H__

#define VERSION           "VERSION 1.5"
#define SIZE_SEED         12
#define NB_DPU            128
#define MAX_NB_DPU_READ   16384                  /* Maximum number of read by DPU by round */
#define MAX_READS_BUFFER  1048576                /* Maximum number of read by round        */
#define MAX_ALIGN         65536

#define COST_SUB         10
#define COST_GAPO        11
#define COST_GAPE        1
#define NB_DIAG          15

#define CODE_SUB 10
#define CODE_DEL 11
#define CODE_INS 12
#define CODE_END 13

typedef struct {
        int size_read;
        int size_neighbour_in_32bits_words;
        int size_neighbour_in_bytes;
        int delta_neighbour_in_bytes;
} reads_info_t;

typedef struct {
        double get_genome;
        double index_genome;
        double get_reads;
        double dispatch_read;
        double map_read;
        double process_read;
        double tot_get_reads;
        double tot_dispatch_read;
        double tot_map_read;
        double tot_process_read;
        double vcf;
} times_ctx_t;

#include <stddef.h>
#include <sys/time.h>
static inline double my_clock(void)
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (1.0e-6 * t.tv_usec + t.tv_sec);
}

#endif /* __UPVC_H__ */
