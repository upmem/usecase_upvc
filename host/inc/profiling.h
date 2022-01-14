#ifndef __PROFILING_H__
#define __PROFILING_H__

#include <time.h>

#define STAT_MAX_SUBSTEPS 10

struct time_stat_t {
    clock_t total_time;
    unsigned int number_calls;
    clock_t substep_total_times[STAT_MAX_SUBSTEPS];
};

struct time_stat_t profiling[14];

#define STAT_DPD                       0
#define STAT_CODE_ALIGNMENT            1
#define STAT_ADD_TO_NON_MAPPED_READ    2
#define STAT_GET_READ_UPDATE_POSITIONS 3
#define STAT_UPDATE_FREQUENCY_TABLE    4
#define STAT_DO_PROCESS_READ           5
#define STAT_PROCESS_READ              6
#define STAT_EXEC_ROUND                7
#define STAT_EXEC_DPUS                 8
#define STAT_THREAD_GET_READS          9
#define STAT_THREAD_DISPATCH           10
#define STAT_THREAD_ACC                11
#define STAT_THREAD_PROCESS            12
#define STAT_DO_MAPPING                13

#define STAT_RECORD_START(FUNCTION)                        \
    clock_t profiling_step_time, profiling_last_step_time; \
    clock_t profiling_start_time = clock();                \
    profiling_last_step_time = profiling_start_time;       \
    profiling[FUNCTION].number_calls++;

#define STAT_RECORD_STEP(FUNCTION, STEP_N)                                                           \
    profiling_step_time = clock();                                                                   \
    profiling[FUNCTION].substep_total_times[STEP_N] += profiling_step_time-profiling_last_step_time; \
    profiling_last_step_time = profiling_step_time;

#define STAT_RECORD_LAST_STEP(FUNCTION, STEP_N)                                      \
    STAT_RECORD_STEP(FUNCTION, STEP_N)                                               \
    profiling[FUNCTION].total_time += profiling_last_step_time-profiling_start_time;


#define PRINT_MICROSECONDS(t)                         \
    if (t<1000000) {                                  \
        printf("%ld.%ldms", t/1000, t%1000);          \
    } else {                                          \
        printf("%ld.%lds", t/1000000, (t/1000)%1000); \
    }

#define PRINT_FUNCTION_STAT(FUNCTION)                                      \
    printf(#FUNCTION ":\n");                                               \
    printf("\tcalled: %u\n", profiling[FUNCTION].number_calls);            \
    printf("\ttotal time:");                                               \
    PRINT_MICROSECONDS(profiling[FUNCTION].total_time)                     \
    printf("\n\tsteps:\n");                                                \
    for (int i=0; i<STAT_MAX_SUBSTEPS; i++)                                \
    {                                                                      \
        if (profiling[FUNCTION].substep_total_times[i] > 0) {              \
            printf("\t\t");                                                \
            PRINT_MICROSECONDS(profiling[FUNCTION].substep_total_times[i]) \
            printf("\n");                                                  \
        }                                                                  \
    }
    

#endif /* __PROFILING_H__ */
