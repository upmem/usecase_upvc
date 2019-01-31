link_directories("${UPMEM_HOME}/usr/lib")

file(GLOB_RECURSE SOURCES src/*.c)
add_executable(upvc ${SOURCES})

target_include_directories(upvc PUBLIC inc/ ../common/inc/)
target_link_libraries(upvc pthread)

add_dpu_dependency(upvc dpu_task)
add_dpu_dependency(upvc dpu_task_opt1)
add_dpu_dependency(upvc dpu_task_opt2)
