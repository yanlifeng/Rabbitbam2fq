cmake_minimum_required(VERSION 3.16)
project(BQC)


set(INC_DIR /Users/ylf9811/ClionProjects/Rabbitbam2fq/htslib-1.11 /usr/local/Cellar/libdeflate/1.6/include)
set(LINK_DIR /Users/ylf9811/ClionProjects/Rabbitbam2fq/HTSLIB/lib /usr/local/Cellar/libdeflate/1.6/lib)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})

set(CMAKE_CXX_STANDARD 14)
add_executable(BQC block_mul.cpp Buffer.cpp Buffer.h BamBlock.cpp BamBlock.h threadconfig.cpp threadconfig.h)
target_link_libraries(BQC libhts.a)
target_link_libraries(BQC libz.a)
target_link_libraries(BQC libdeflate.a)