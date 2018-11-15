INCLUDE_DIRECTORIES(inc)
INCLUDE_DIRECTORIES(../common/inc)

file(GLOB_RECURSE includes inc/*.h)
message("includes files: ${includes}")

file(GLOB_RECURSE common_includes ../common/inc/*.h)
message("common includes files: ${common_includes}")

set(SOURCES_OPT0 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd.c src/task.c)
dpu_generate_configuration(dpu_task_config.json build_configuration.script)
dpu_build(dpu_task dpu_task.bin dpu_task_config.json -flto "${SOURCES_OPT0}" "${common_includes};${includes}")

set(SOURCES_OPT1 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt1.c src/task.c)
dpu_generate_configuration(dpu_task_config_opt1.json build_configuration_opt1.script)
dpu_build(dpu_task_opt1 dpu_task_opt1.bin dpu_task_config_opt1.json "" "${SOURCES_OPT1}" "${common_includes};${includes}")

set(SOURCES_OPT2 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd_init_opt.c src/odpd_opt.s src/task.c)
dpu_generate_configuration(dpu_task_config_opt2.json build_configuration_opt2.script)
dpu_build(dpu_task_opt2 dpu_task_opt2.bin dpu_task_config_opt2.json "" "${SOURCES_OPT2}" "${common_includes};${includes}")
