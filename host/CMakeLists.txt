cmake_minimum_required(VERSION 2.8.12)
project(upvc)

if (NOT DEFINED UPMEM_HOME)
        if ( "$ENV{UPMEM_HOME}" STREQUAL "")
                set(UPMEM_HOME "/usr")
        else ()
                set(UPMEM_HOME $ENV{UPMEM_HOME})
        endif ()
endif ()

include(${UPMEM_HOME}/share/upmem/cmake/include/host/DpuHost.cmake)

set(CMAKE_C_FLAGS "--std=c99 -O3 -Wall -Wextra -Werror -g3")
link_directories("${DPU_HOST_LINK_DIRECTORIES}")

file(GLOB_RECURSE SOURCES src/*.c)
add_executable(upvc ${SOURCES})

target_include_directories(upvc PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/)
target_link_libraries(upvc ${DPU_HOST_LIBRARIES} pthread)

set(DPU_PROJECT_RELATIVE_PATH ../dpu)

if(NOT DPU_TARGET)
        set(DPU_TARGET "fsim")
        if(NOT NB_DPU)
                set(NB_DPU "2")
        endif()
endif()

set(NB_DPU_MARK)
if (NB_DPU)
        set(NB_DPU_MARK "-n")
endif()

add_custom_target(run ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map ${NB_DPU_MARK} ${NB_DPU} -t ${DPU_TARGET} -b ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}/dpu_task
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
        CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/share/upmem/cmake/dpu.cmake -DUPMEM_HOME=${UPMEM_HOME}
        BUILD_ALWAYS TRUE
        INSTALL_COMMAND ""
)

add_dependencies(upvc dpu_task)
