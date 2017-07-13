cmake_minimum_required(VERSION 3.7)
project(GSM)

set(CMAKE_CXX_STANDARD 11)

include_directories(C:/Eigen/Eigen)
include_directories(X:/Eigen/Eigen)
include_directories(${CMAKE_SOURCE_DIR}/src)

set(SOURCE_FILES main.cpp src/vector.h src/GSM_Solver.h src/GSM_Solver.cpp src/GA_GSM_Solver.h)
add_executable(GSM ${SOURCE_FILES})