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
target_link_libraries(upvc dpu pthread)

set(DPU_PROJECT_RELATIVE_PATH ../dpu)

if(NOT NB_DPU)
        set(NB_DPU "2")
endif()

if(NOT DPU_TARGET)
        set(DPU_TARGET "hsim")
endif()

add_custom_target(run ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map -n ${NB_DPU} -t ${DPU_TARGET} -b ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}/dpu_task
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

add_custom_target(check diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta.ref
        COMMAND diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta.ref
        COMMAND diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta.ref
        COMMAND diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta.ref
        COMMAND diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta.ref
        COMMAND diff ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta.ref
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

include(ExternalProject)
ExternalProject_Add(
        dpu_task
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/usr/share/upmem/cmake/dpu.cmake -DUPMEM_HOME=${UPMEM_HOME}
        INSTALL_COMMAND ""
)

add_dependencies(upvc dpu_task)
