cmake_minimum_required(VERSION 3.0)
project(upvc)

set(CMAKE_C_FLAGS "--std=c99 -O3 -Wall -Wextra -Werror -g3")

set(DPU_DIRECTORY dpu)
set(CMAKE_UTILS_DIRECTORY cmake_utils)
include(cmake_utils/CMakeDpu.cmake)

add_subdirectory(host)
add_subdirectory(dpu)
