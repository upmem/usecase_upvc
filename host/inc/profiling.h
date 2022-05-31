#ifndef __PROFILING_H__
#define __PROFILING_H__

#include <time.h>

#define PROFILE_STAT false
#define STAT_MAX_SUBSTEPS 10

struct time_stat_t {
    clock_t total_time;
    unsigned int number_calls;
    clock_t substep_total_times[STAT_MAX_SUBSTEPS];
};

struct time_stat_t profiling[25];

#define STAT_SUB_ONLY_PATH             0
#define STAT_DPD                       1
#define STAT_CODE_ALIGNMENT            2
#define STAT_ADD_TO_NON_MAPPED_READ    3
#define STAT_GET_READ_UPDATE_POSITIONS 4
#define STAT_UPDATE_FREQUENCY_TABLE    5
#define STAT_DO_PROCESS_READ           6
#define STAT_PROCESS_READ              7
#define STAT_EXEC_ROUND                8
#define STAT_EXEC_DPUS                 9
#define STAT_THREAD_GET_READS          10
#define STAT_THREAD_DISPATCH           11
#define STAT_THREAD_ACC                12
#define STAT_THREAD_PROCESS            13
#define STAT_DO_MAPPING                14
#define STAT_ADD_CODEPENDENCE_INFO     15
#define STAT_GET_NEW_CODEPENDENCE_INFO 16
#define STAT_ALLOCATE_NEW_CHUNK        17

#define PRINT_MICROSECONDS(t)                         \
    if (t<CLOCKS_PER_SEC) {                           \
        printf("%fms", (float)t/CLOCKS_PER_SEC*1000); \
    } else {                                          \
        printf("%fs", (float)t/CLOCKS_PER_SEC);       \
    }


#if PROFILE_STAT

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


#define PRINT_FUNCTION_STAT(FUNCTION)                                      \
    printf(#FUNCTION ":\n");                                               \
    printf("\tcalled: %u\n", profiling[FUNCTION].number_calls);            \
    printf("\ttotal time:");                                               \
    PRINT_MICROSECONDS(profiling[FUNCTION].total_time)                     \
    printf("\n\tsteps:\n");                                                \
    for (int i=0; i<STAT_MAX_SUBSTEPS; i++)                                \
    {                                                                      \
        if (profiling[FUNCTION].substep_total_times[i] > 0) {              \
            printf("\t\t%d:\t", i);                                        \
            PRINT_MICROSECONDS(profiling[FUNCTION].substep_total_times[i]) \
            printf("\n");                                                  \
        }                                                                  \
    }

#define PRINT_ALL_FUNCTION_STAT()                                       \
    printf("\nprofiling:\n\n");                                         \
    PRINT_FUNCTION_STAT(STAT_DO_MAPPING);                               \
    PRINT_FUNCTION_STAT(STAT_THREAD_PROCESS);                           \
    PRINT_FUNCTION_STAT(STAT_THREAD_ACC);                               \
    PRINT_FUNCTION_STAT(STAT_THREAD_DISPATCH);                          \
    PRINT_FUNCTION_STAT(STAT_THREAD_GET_READS);                         \
    PRINT_FUNCTION_STAT(STAT_EXEC_DPUS);                                \
    PRINT_FUNCTION_STAT(STAT_EXEC_ROUND);                               \
    PRINT_FUNCTION_STAT(STAT_PROCESS_READ);                             \
    PRINT_FUNCTION_STAT(STAT_DO_PROCESS_READ);                          \
    PRINT_FUNCTION_STAT(STAT_ADD_TO_NON_MAPPED_READ);                   \
    PRINT_FUNCTION_STAT(STAT_UPDATE_FREQUENCY_TABLE);                   \
    PRINT_FUNCTION_STAT(STAT_GET_READ_UPDATE_POSITIONS);                \
    PRINT_FUNCTION_STAT(STAT_CODE_ALIGNMENT);                           \
    PRINT_FUNCTION_STAT(STAT_SUB_ONLY_PATH);                            \
    PRINT_FUNCTION_STAT(STAT_DPD);                                      \
    PRINT_FUNCTION_STAT(STAT_ADD_CODEPENDENCE_INFO);                    \
    PRINT_FUNCTION_STAT(STAT_GET_NEW_CODEPENDENCE_INFO);                \
    PRINT_FUNCTION_STAT(STAT_ALLOCATE_NEW_CHUNK);

#else

#define STAT_RECORD_START(FUNCTION)
#define STAT_RECORD_STEP(FUNCTION, STEP_N)
#define STAT_RECORD_LAST_STEP(FUNCTION, STEP_N)
#define PRINT_FUNCTION_STAT(FUNCTION)
#define PRINT_ALL_FUNCTION_STAT()

#endif
    

#endif /* __PROFILING_H__ */
