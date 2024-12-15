cmake_minimum_required(VERSION 3.12)
project(mpna-matrix-storage C)

set(CMAKE_C_STANDARD 11)

add_executable(store main.c)
target_compile_options(store PRIVATE "-O3" "-Wall" "-Wno-incompatible-pointer-types")