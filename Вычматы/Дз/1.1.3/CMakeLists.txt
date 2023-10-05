cmake_minimum_required(VERSION 3.22)
project(1_1_2)

set(CMAKE_CXX_STANDARD 17)

add_executable(1_1_2 newton.cpp lesha.cpp)
add_executable(second spline.cpp)
