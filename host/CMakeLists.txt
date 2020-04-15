cmake_minimum_required(VERSION 2.8.12)
project(upvc)

set(CMAKE_C_FLAGS "--std=gnu99 -O3 -Wall -Wextra -Werror -g3 -fopenmp")

file(GLOB_RECURSE SOURCES src/*.c)
add_executable(upvc ${SOURCES})

target_include_directories(upvc PUBLIC inc/ ../common/inc/)
target_link_libraries(upvc pthread)

set(NB_DPU_MARK)
if (NB_DPU)
        set(NB_DPU_MARK "-n")
endif()

add_custom_target(run ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map ${NB_DPU_MARK} ${NB_DPU}
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

add_custom_target(check cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta.ref | sort > ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta.ref.sort
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE1.fasta.ref.sort -
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE1.fasta.ref.sort -
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE1.fasta.ref.sort -
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_1_PE2.fasta.ref.sort -
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_2_PE2.fasta.ref.sort -
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta | sort | diff -y ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_3_PE2.fasta.ref.sort -
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
