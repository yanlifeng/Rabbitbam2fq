cmake_minimum_required(VERSION 3.16)
project(Rabbitbam2fq)

#set(INC_DIR /Users/ylf9811/ClionProjects/Rabbitbam2fq/HTSLIB/include)
set(INC_DIR /Users/ylf9811/ClionProjects/Rabbitbam2fq/htslib-1.11)
set(LINK_DIR /Users/ylf9811/ClionProjects/Rabbitbam2fq/HTSLIB/lib)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(Rabbitbam2fq test2bt.c)
target_link_libraries(Rabbitbam2fq libhts.a)
