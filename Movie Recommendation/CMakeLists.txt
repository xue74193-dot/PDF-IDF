cmake_minimum_required(VERSION 3.17)
project(Mov)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(Mov haba_non.h main.cpp utility_functions.h RandomMultiGreedy.h ParallelDistortedFilter.h Boost.h Boost_plus.h RandomSelection.h IterativeDistortedFilter.h)
