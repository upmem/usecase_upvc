cmake_minimum_required(VERSION 3.0)
project(upvc)

if (NOT DEFINED UPMEM_HOME)
    set(UPMEM_HOME $ENV{UPMEM_HOME})
endif ()

set(CMAKE_C_FLAGS "--std=c99 -O3 -Wall -Wextra -Werror -g3")
link_directories("${UPMEM_HOME}/usr/lib")

file(GLOB_RECURSE SOURCES src/*.c)
add_executable(upvc ${SOURCES})

target_include_directories(upvc PUBLIC ${UPMEM_HOME}/usr/share/upmem/include/host inc/ ../common/inc/)
target_link_libraries(upvc pthread dpu dpucni)

set(DPU_PROJECT_RELATIVE_PATH ../dpu)
include(ExternalProject)
ExternalProject_Add(
        dpu_task
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/usr/share/upmem/cmake/dpu.cmake -DCMAKE_BUILD_TYPE=MinSizeRel -DUPMEM_HOME=${UPMEM_HOME}
        INSTALL_COMMAND ""
)

add_dependencies(upvc dpu_task)
