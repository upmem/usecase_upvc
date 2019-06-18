/**
 * @Copyright (c) 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __STATS_H__
#define __STATS_H__

#ifndef STATS_ON

#define STATS_ATTRIBUTE __attribute__((unused))
#define STATS_INCR_NB_REQS(stats)
#define STATS_INCR_NB_NODP_CALLS(stats)
#define STATS_INCR_NB_ODPD_CALLS(stats)
#define STATS_INCR_NB_RESULTS(stats, val)
#define STATS_INCR_LOAD(stats_ptr, val)
#define STATS_INCR_LOAD_DATA(stats_ptr, val)
#define STATS_INCR_STORE(stats_ptr, val)
#define STATS_INCR_STORE_RESULT(stats_ptr, val)
#define STATS_GET_START_TIME(start, acc, end)
#define STATS_GET_END_TIME(end, acc)
#define STATS_STORE_NODP_TIME(stats_ptr, val)
#define STATS_STORE_ODPD_TIME(stats_ptr, val)
#define STATS_TIME_VAR(start, end, acc)

#else /* STATS_ON */

#define STATS_ATTRIBUTE
#define STATS_INCR_NB_REQS(stats)                                                                                                \
    do {                                                                                                                         \
        (stats).nb_reqs++;                                                                                                       \
    } while (0)
#define STATS_INCR_NB_NODP_CALLS(stats)                                                                                          \
    do {                                                                                                                         \
        (stats).nb_nodp_calls++;                                                                                                 \
    } while (0)
#define STATS_INCR_NB_ODPD_CALLS(stats)                                                                                          \
    do {                                                                                                                         \
        (stats).nb_odpd_calls++;                                                                                                 \
    } while (0)
#define STATS_INCR_NB_RESULTS(stats, val)                                                                                        \
    do {                                                                                                                         \
        (stats).nb_results += (val);                                                                                             \
    } while (0)
#define STATS_INCR_LOAD(stats_ptr, val)                                                                                          \
    do {                                                                                                                         \
        (stats_ptr)->mram_load += (val);                                                                                         \
    } while (0)
#define STATS_INCR_LOAD_DATA(stats_ptr, val)                                                                                     \
    do {                                                                                                                         \
        (stats_ptr)->mram_data_load += (val);                                                                                    \
    } while (0)
#define STATS_INCR_STORE(stats_ptr, val)                                                                                         \
    do {                                                                                                                         \
        (stats_ptr)->mram_store += (val);                                                                                        \
    } while (0)
#define STATS_INCR_STORE_RESULT(stats_ptr, val)                                                                                  \
    do {                                                                                                                         \
        (stats_ptr)->mram_result_store += (val);                                                                                 \
    } while (0)
#define STATS_GET_START_TIME(start, acc, end)                                                                                    \
    do {                                                                                                                         \
        start = perfcounter_get();                                                                                               \
        acc = 0ULL;                                                                                                              \
        end = start;                                                                                                             \
    } while (0)
#define STATS_GET_END_TIME(end, acc)                                                                                             \
    do {                                                                                                                         \
        get_time_and_accumulate(&acc, &end);                                                                                     \
    } while (0)
#define STATS_STORE_NODP_TIME(stats_ptr, val)                                                                                    \
    do {                                                                                                                         \
        (stats_ptr)->nodp_time += (val);                                                                                         \
    } while (0)
#define STATS_STORE_ODPD_TIME(stats_ptr, val)                                                                                    \
    do {                                                                                                                         \
        (stats_ptr)->odpd_time += (val);                                                                                         \
    } while (0)
#define STATS_TIME_VAR(start, end, acc) perfcounter_t start, end, acc;

#endif /* STATS_ON */
#endif /* __STATS_H__ */
