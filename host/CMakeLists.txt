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

set(DPU_PROJECT_RELATIVE_PATH ../dpu)
set(DPU_BINARY_NAME dpu_task)
set(NR_TASKLETS 16)
set(CMAKE_C_FLAGS "--std=gnu99 -O3 -Wall -Wextra -Werror -g3 -DNR_TASKLETS=${NR_TASKLETS} -fopenmp -DDPU_BINARY=\\\"${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}/${DPU_BINARY_NAME}\\\"")
link_directories("${DPU_HOST_LINK_DIRECTORIES}")

file(GLOB_RECURSE SOURCES src/*.c)
add_executable(upvc ${SOURCES})

target_include_directories(upvc PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/)
target_link_libraries(upvc ${DPU_HOST_LIBRARIES} pthread)

set(NB_DPU_MARK)
if (NB_DPU)
        set(NB_DPU_MARK "-n")
endif()

if (NOT DEFINED INDEX_WITH_DPUS)
        set(INDEX_WITH_DPUS "")
endif()

add_custom_target(mrams ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g index -n 128 ${INDEX_WITH_DPUS}
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(run ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map ${NB_DPU_MARK} ${NB_DPU} -f
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(run_simu ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map -f -s
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

add_custom_target(check cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc.vcf | grep -E "^chr22" | diff -y --suppress-common-lines ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc_ref.vcf -
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

include(ExternalProject)
ExternalProject_Add(
        ${DPU_BINARY_NAME}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/share/upmem/cmake/dpu.cmake -DUPMEM_HOME=${UPMEM_HOME} -DNR_TASKLETS=${NR_TASKLETS}
        BUILD_ALWAYS TRUE
        INSTALL_COMMAND ""
)

add_dependencies(upvc dpu_task)
