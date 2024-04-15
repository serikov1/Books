cmake_minimum_required(VERSION 3.26)
project(2_2_2)

set(CMAKE_CXX_STANDARD 17)

add_executable(task main.cpp)
add_subdirectory(fmt-10.2.0)
find_package(fmt)
target_link_libraries(task fmt::fmt)