cmake_minimum_required(VERSION 3.6)
project(kmeans)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Werror")

add_executable(kmeans main.cpp)