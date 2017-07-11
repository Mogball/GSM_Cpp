cmake_minimum_required(VERSION 3.7)
project(GSM)

set(CMAKE_CXX_STANDARD 11)

include_directories(C:/Eigen/Eigen)

set(SOURCE_FILES main.cpp)
add_executable(GSM ${SOURCE_FILES})