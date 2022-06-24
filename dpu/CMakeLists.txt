cmake_minimum_required(VERSION 2.8.12)
project(dpu C ASM)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)

set(CMAKE_C_FLAGS "-O2 -g -fstack-size-section -DNR_TASKLETS=${NR_TASKLETS} -DSTACK_SIZE_DEFAULT=192")

INCLUDE_DIRECTORIES(inc)
INCLUDE_DIRECTORIES(../common/inc/)

set(SOURCES_OPT2 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt.S src/task.c src/nodp.S)

add_executable(dpu_task ${SOURCES_OPT2})
