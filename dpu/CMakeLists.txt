INCLUDE_DIRECTORIES(inc)

file(GLOB_RECURSE includes inc/*.h)
message("includes files: ${includes}")

set(SOURCES_OPT0 src/odpd.c src/task.c)
dpu_generate_configuration(dpu_task_config.json build_configuration.script)
dpu_build(dpu_task dpu_task.bin dpu_task_config.json -flto "${SOURCES_OPT0}" "${includes}")

set(SOURCES_OPT1 src/odpd_opt1.c src/task.c)
dpu_generate_configuration(dpu_task_config_opt1.json build_configuration_opt1.script)
dpu_build(dpu_task_opt1 dpu_task_opt1.bin dpu_task_config_opt1.json "" "${SOURCES_OPT1}" "${includes}")

set(SOURCES_OPT2 src/odpd_opt.s src/odpd_opt1.c src/task.c)
dpu_generate_configuration(dpu_task_config_opt2.json build_configuration_opt2.script)
dpu_build(dpu_task_opt2 dpu_task_opt2.bin dpu_task_config_opt2.json "-DOPT_S" "${SOURCES_OPT2}" "${includes}")
