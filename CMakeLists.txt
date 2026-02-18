cmake_minimum_required(VERSION 4.1)
project(Computing C)

set(CMAKE_C_STANDARD 23)

add_executable(Computing main.c
        ode_utils.h
        ode_utils.c)
