cmake_minimum_required(VERSION 3.0)
project(tests_tools)

set(CMAKE_C_FLAGS "--std=c99 -O3 -Wall -Werror -Wextra -g3")

add_executable(upvclog2graph upvclog2graph.c)
add_executable(csv2svg csv2svg.c)
add_executable(extract_res extract_res.c)
add_executable(extract_trace extract_trace.c)
add_executable(seed_txt_2_seed_raw seed_txt_2_seed_raw.c)
