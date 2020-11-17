cmake_minimum_required(VERSION 3.17)
project(kmeans C)

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi -Wall -Wextra -Werror -pedantic-errors")
add_executable(kmeans kmeans.c)