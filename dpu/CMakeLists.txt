cmake_minimum_required(VERSION 3.0)
project(dpu C ASM)

INCLUDE_DIRECTORIES(inc)
INCLUDE_DIRECTORIES(../common/inc/)

include(DpuTools)

# set(SOURCES_OPT0 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd.c src/task.c)
# dpu_build_bootstrap(dpu_task_config.S SCRIPT build_configuration.script)

# set(SOURCES_OPT1 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt1.c src/task.c)
# dpu_build_bootstrap(dpu_task_config.S SCRIPT build_configuration_opt1.script)

set(SOURCES_OPT2 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt.S src/task.c)
dpu_build_bootstrap(dpu_task_config.S SCRIPT build_configuration_opt2.script)

add_executable(dpu_task dpu_task_config.S ${SOURCES_OPT2})
