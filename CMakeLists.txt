cmake_minimum_required(VERSION 3.27)
project(ComputationalMethod)

set(CMAKE_CXX_STANDARD 14)

add_executable(ComputationalMethod
        hw6/src/main.cpp
        hw6/src/FastFourier.cpp
        hw6/src/FastFourier.h)
