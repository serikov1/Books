cmake_minimum_required(VERSION 3.26)
project(1_2_1)

set(CMAKE_CXX_STANDARD 17)

add_executable(1_2_1 main.cpp)
add_executable(second second.cpp)
add_subdirectory(fmt-master)
find_package(fmt)
target_link_libraries(1_2_1 fmt::fmt)
target_link_libraries(second fmt::fmt)