cmake_minimum_required(VERSION 2.8.12)
project(tests_tools)

set(CMAKE_C_FLAGS "--std=c99 -O3 -Wall -Werror -Wextra -g3")

add_executable(upvclog2graph upvclog2graph.c)
add_executable(csv2svg csv2svg.c)
add_executable(extract_res extract_res.c)
