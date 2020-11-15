cmake_minimum_required(VERSION 3.17)
project(K_MEANS C)

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Werror -pedantic-errors")
add_executable(K_MEANS main.c)