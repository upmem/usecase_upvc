#define V_QUIET 0
#define V_FATAL 1
#define V_ERROR 2
#define V_WARN 3
#define V_INFO 4
#define V_DEBUG 5
#define V_TRACE 6

#define VERBOSE V_TRACE
#define VERBOSE_COLORS true
#define VERBOSE_LOG_LEVEL true

#if VERBOSE_COLORS
#define VERBOSE_COLOR_START_FATAL "\033[41m"
#define VERBOSE_COLOR_START_ERROR "\033[31m"
#define VERBOSE_COLOR_START_WARN  "\033[33m"
#define VERBOSE_COLOR_START_INFO  "\033[32m"
#define VERBOSE_COLOR_START_DEBUG "\033[34m"
#define VERBOSE_COLOR_START_TRACE "\033[36m"
#define VERBOSE_COLOR_END         "\033[0m"
#else
#define VERBOSE_COLOR_START_FATAL
#define VERBOSE_COLOR_START_ERROR
#define VERBOSE_COLOR_START_WARN
#define VERBOSE_COLOR_START_INFO
#define VERBOSE_COLOR_START_DEBUG
#define VERBOSE_COLOR_START_TRACE
#define VERBOSE_COLOR_END
#endif

#if VERBOSE_LOG_LEVEL
#define VERBOSE_PRINT_PREFIX(level) VERBOSE_COLOR_START_##level #level "\t" VERBOSE_COLOR_END
#else
#define VERBOSE_PRINT_PREFIX(level)
#endif


#if VERBOSE>=V_TRACE
#define LOG_TRACE(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(TRACE) __VA_ARGS__)
#else
#define LOG_TRACE(...)
#endif

#if VERBOSE>=V_DEBUG
#define LOG_DEBUG(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(DEBUG) __VA_ARGS__)
#else
#define LOG_DEBUG(...)
#endif

#if VERBOSE>=V_INFO
#define LOG_INFO(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(INFO) __VA_ARGS__)
#else
#define LOG_INFO(...)
#endif

#if VERBOSE>=V_WARN
#define LOG_WARN(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(WARN) __VA_ARGS__)
#else
#define LOG_WARN(...)
#endif

#if VERBOSE>=V_ERROR
#define LOG_ERROR(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(ERROR) __VA_ARGS__)
#else
#define LOG_ERROR(...)
#endif

#if VERBOSE>=V_FATAL
#define LOG_FATAL(...) fprintf(stderr, VERBOSE_PRINT_PREFIX(FATAL) __VA_ARGS__)
#else
#define LOG_FATAL(...)
#endif
