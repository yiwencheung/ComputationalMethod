cmake_minimum_required(VERSION 3.27)
project(ComputationalMethod)

set(CMAKE_CXX_STANDARD 14)

add_executable(ComputationalMethod
        hw3/src/main.cpp
        hw3/src/InversePowerMethod.cpp
        hw3/src/InversePowerMethod.h)
