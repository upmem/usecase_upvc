/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_DEBUG_H__
#define __INTEGRATION_DEBUG_H__

#include "common.h"

/* Define DEBUG to get plain-text traces */
/* #define DEBUG */

/* To print out every single request (pretty slow) */
/* #define DEBUG_REQUESTS */

/* To print the preliminary summary of MRAM mapping */
/* #define DEBUG_MRAM_INFO */

/* To print the intermediate computation results */
/* #define DEBUG_RESULTS */

/* To print the different calls and occurrences in the main loop */
/* #define DEBUG_STATS */

/* To print the load/store MRAM performed by each tasklet */
/* #define DEBUG_STATS_MRAM */

/* Define NB_RUNNING_TASKLETS to define the effective number of operating tasklets (1 to 16) */
/* #define NB_RUNNING_TASKLETS 1 */

/* Define ASSERT_DMA_ADDR_ALIGNMENT to verify that the base addresses of a DMA transfer are */
/* properly aligned on 64 bits */
/* #define ASSERT_DMA_ADDR_ALIGNMENT */

/* To generate reference patterns for internal developments */
/* #define DEBUG_PROCESS */




#ifdef DEBUG
#include <ktrace.h>
#define printf ktrace
#else /* DEBUG */
#define printf(...)
#endif /* DEBUG */




#ifndef NB_RUNNING_TASKLETS
#define NB_RUNNING_TASKLETS NB_TASKLET_PER_DPU
#endif /* NB_RUNNING_TASKLETS */




#ifdef DEBUG_MRAM_INFO
#define DEBUG_MRAM_INFO_PRINT(mram_info)                                \
        do {                                                            \
                printf("MRAM total_nbr_size = %u\n", (mram_info)->total_nbr_size); \
                printf("MRAM nb neigbors    = %u\n", (mram_info)->nb_nbr); \
                printf("MRAM nbr len        = %u\n", (mram_info)->nbr_len); \
                printf("MRAM magic          = 0x%x\n", (mram_info)->magic); \
        } while(0)
#else /* DEBUG_MRAM_INFO */
#define DEBUG_MRAM_INFO_PRINT(mram_info)
#endif /* DEBUG_MRAM_INFO */




#ifdef DEBUG_RESULTS

#define DEBUG_RESULTS_VAR int debug_nr = -1
#define DEBUG_RESULTS_INCR (debug_nr++)
#define DEBUG_RESULTS_PRINT(thread, idx, request_num, request_offset, seed, seq, score) \
        do {                                                            \
                printf("[%u] nr=%u ix=%u num=%u ", (thread), debug_nr, (idx), (request_num)); \
                printf("offset=%u seed=%u seq=%u score=%u\n",           \
                       (request_offset),                                \
                       (seed),                                          \
                       (seq),                                           \
                       (score));                                        \
        } while(0)

#else /* DEBUG_RESULTS */

#define DEBUG_RESULTS_VAR
#define DEBUG_RESULTS_INCR
#define DEBUG_RESULTS_PRINT(thread, idx, request_num, request_offset, seed, seq, score)

#endif /* DEBUG_RESULTS */




#ifdef ASSERT_DMA_ADDR_ALIGNMENT

#include <defs.h>
#define ASSERT_DMA_ADDR_ALIGNMENT_MRAM_SIZE (64 << 20)
#define ASSERT_DMA_ADDR_ALIGNMENT_WRAM_SIZE (128 << 10)
#define ASSERT_DMA_ADDR(mram_addr, wram_addr, len)                      \
        do {                                                            \
                if (((mram_addr) + (len)) > ASSERT_DMA_ADDR_ALIGNMENT_MRAM_SIZE) { \
                        printf("MRAM buffer [%x,%x[ does not fit in MRAM!\n", (mram_addr), (mram_addr) + (len)); \
                        halt();                                         \
                }                                                       \
                if (((unsigned int) (wram_addr) + (len)) > ASSERT_DMA_ADDR_ALIGNMENT_WRAM_SIZE) { \
                        printf("WRAM buffer [%x,%x[ does not fit in WRAM!\n", (wram_addr), (wram_addr) + (len)); \
                        halt();                                         \
                }                                                       \
                if (((mram_addr) & 7) != 0) {                           \
                        printf("MRAM address %x is not aligned on longs!\n", (mram_addr)); \
                        halt();                                         \
                }                                                       \
                if ((((unsigned int) (wram_addr)) & 7) != 0) {          \
                        printf("WRAM address %x is not aligned on longs!\n", (wram_addr)); \
                        halt();                                         \
                }                                                       \
        } while(0)
#define ASSERT_DMA_LEN(len)                                             \
        do {                                                            \
                if (((len) & 7) != 0) {                                 \
                        printf("DMA transfer len %u is not a multiple of longs!\n", (len)); \
                        halt();                                         \
                }                                                       \
        } while(0)

#else

#define ASSERT_DMA_ADDR(x, y, l)
#define ASSERT_DMA_LEN(x)

#endif /* ASSERT_DMA_ALIGNMENT */




#ifdef DEBUG_STATS

 /**
 * @var nb_reads       Number of reads processed by this.
 * @var nb_refs        Number of references processed by this.
 * @var nb_odpd_calls  Number of calls to ODPD.
 */
typedef struct stats {
        unsigned int nb_reads;
        unsigned int nb_refs;
        unsigned int nb_odpd_calls;
} *stats_t;

#define DEBUG_STATS_INC_NB_READS do { stats.nb_reads++;} while(0)
#define DEBUG_STATS_INC_NB_REFS do { stats.nb_refs++;} while(0)
#define DEBUG_STATS_INC_NB_ODPD_CALLS  do {stats.nb_odpd_calls++;} while(0)
#define DEBUG_STATS_PRINT                                               \
        do {                                                            \
                printf("STATS nb_reads=%u nb_refs=%u nb_odpd_calls=%u\n", \
                       stats.nb_reads,                                  \
                       stats.nb_refs,                                   \
                       stats.nb_odpd_calls);                            \
        } while(0)
#define DEBUG_STATS_CLEAR                       \
        do {                                    \
                stats.nb_reads = 0;             \
                stats.nb_refs = 0;              \
                stats.nb_odpd_calls = 0;        \
        } while(0)
#define DEBUG_STATS_STRUCT struct stats stats

#else

#define DEBUG_STATS_INC_NB_READS
#define DEBUG_STATS_INC_NB_REFS
#define DEBUG_STATS_INC_NB_ODPD_CALLS
#define DEBUG_STATS_PRINT
#define DEBUG_STATS_CLEAR
#define DEBUG_STATS_STRUCT

#endif /* DEBUG_STATS */




#ifdef DEBUG_REQUESTS

#define DEBUG_REQUESTS_VAR \
        uint8_t *debug_process_syms = mem_alloc_dma(NB_BYTES_TO_SYMS(ALIGN_DPU(mram_info.nbr_len)));
#define DEBUG_REQUESTS_PRINT_POOL(request_pool)                         \
        do {                                                            \
                printf("R> nb_reads = %u\n", (request_pool).nb_reads);  \
                printf("R> rdidx    = %u\n", (request_pool).rdidx);     \
                printf("R> cur_read = %x\n", (uint32_t) ((request_pool).cur_read)); \
        } while(0)
#define DEBUG_REQUESTS_DECODE(bytes, nbr_len)                           \
        do {                                                            \
                for (int i = 0; i < NB_BYTES_TO_SYMS(nbr_len); i += 4) { \
                        uint8_t this_byte = (bytes)[i >> 2];            \
                        debug_process_syms[i] = (uint8_t) (this_byte & 3); \
                        debug_process_syms[i + 1] = (uint8_t) ((this_byte >> 2) & 3); \
                        debug_process_syms[i + 2] = (uint8_t) ((this_byte >> 4) & 3); \
                        debug_process_syms[i + 3] = (uint8_t) ((this_byte >> 6) & 3); \
                }                                                       \
        } while (0)
#define DEBUG_REQUESTS_PRINT_SYMS(nbr_len)                              \
        do {                                                            \
                char as_char;                                           \
                for (int i = 0; i < NB_BYTES_TO_SYMS(nbr_len); i++) {   \
                        switch (debug_process_syms[i]) {                \
                        case 0:                                         \
                                as_char = 'A';                          \
                                break;                                  \
                        case 1:                                         \
                                as_char = 'C';                          \
                                break;                                  \
                        case 2 :                                        \
                                as_char = 'T';                          \
                                break;                                  \
                        default:                                        \
                                as_char = 'G';                          \
                        }                                               \
                        printf("%c", as_char);                          \
                }                                                       \
                printf("\n");                                           \
        } while (0)
#define DEBUG_REQUESTS_PRINT(request, nbr_len)                          \
        do {                                                            \
                printf("request: offset %u count %u num %u\n", (request)->offset, (request)->count, (request)->num); \
                DEBUG_REQUESTS_DECODE((uint8_t *)(((uint8_t*)(request))+sizeof(dpu_request_t)), nbr_len); \
                DEBUG_REQUESTS_PRINT_SYMS(nbr_len);                     \
        } while(0)
#define DEBUG_REQUESTS_PRINT_REF(nbr, nbr_len)                          \
        do {                                                            \
                printf("ref\n");                                        \
                DEBUG_REQUESTS_DECODE(nbr, nbr_len);                    \
                DEBUG_REQUESTS_PRINT_SYMS(nbr_len);                     \
        } while(0)

#else

#define DEBUG_REQUESTS_VAR
#define DEBUG_REQUESTS_PRINT_POOL(request_pool)
#define DEBUG_REQUESTS_PRINT_REF(request, nbr_len)
#define DEBUG_REQUESTS_PRINT(request, nbr_len)

#endif /* DEBUG_REQUESTS */




#ifdef DEBUG_PROCESS

#define DEBUG_PROCESS_PROFILE(len) do { printf("NBRLEN:%u\n", len);} while(0)
#define DEBUG_PROCESS_PATTERN(name, nb_bytes, req)                      \
        do {                                                            \
                unsigned int each_byte;                                 \
                printf("%s:", name);                                    \
                for (each_byte = 0; each_byte < nb_bytes; each_byte++) { \
                        printf(" %02x", req[each_byte]);                \
                }                                                       \
                printf("\n");                                           \
        } while(0)
#define DEBUG_PROCESS_SCORES(nb_bytes, req, ref, mini, score_nodp, score_odpd, lockit) \
        do {                                                            \
                mutex_lock(lockit);                                     \
                DEBUG_PROCESS_PATTERN("REQ", nb_bytes, req);            \
                DEBUG_PROCESS_PATTERN("REF", nb_bytes, ref);            \
                printf("mini %u\n", mini);                              \
                printf("score nodp %u\n", score_nodp);                  \
                if (score_odpd != -1) {                                 \
                        printf("score odpd %u\n", score_odpd);          \
                }                                                       \
                mutex_unlock(lockit);                                   \
        } while(0)

#else

#define DEBUG_PROCESS_PROFILE(n)
#define DEBUG_PROCESS_SCORES(n, req, ref, mini, snodp, sodpd, lockit)

#endif /* TRACE_PATTERNS_AND_SCORES */




#ifdef DEBUG_STATS_MRAM
#define DEBUG_STATS_MRAM_PRINT(stats_load, stats_store, lockit)         \
        do {                                                            \
                mutex_lock(lockit);                                     \
                printf("mem stats: reads=%d writes=%d\n", (stats_load), (stats_store)); \
                mutex_unlock(lockit);                                   \
        } while(0)
#else
#define DEBUG_STATS_MRAM_PRINT(stats_load, stats_store, lockit)
#endif /* DEBUG_STATS_MRAM */

#endif /* __INTEGRATION_DEBUG_H__ */
