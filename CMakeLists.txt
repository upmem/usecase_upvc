cmake_minimum_required(VERSION 3.6)
project(upvc)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -Werror -g3")

set(DPU_DIRECTORY dpu)
set(CMAKE_UTILS_DIRECTORY cmake_utils)
include(cmake_utils/CMakeDpu.cmake)

add_subdirectory(host)
add_subdirectory(dpu)
