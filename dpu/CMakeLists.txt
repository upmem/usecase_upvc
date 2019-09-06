cmake_minimum_required(VERSION 2.8.12)
project(dpu C ASM)

set(CMAKE_C_FLAGS -O2)

INCLUDE_DIRECTORIES(inc)
INCLUDE_DIRECTORIES(../common/inc/)

set(SOURCES_OPT2 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt.S src/task.c)

add_executable(dpu_task ${SOURCES_OPT2})
