cmake_minimum_required(VERSION 3.20)
project(rtfmm_3d VERSION 1.0 LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 20)

message("configure all")

add_subdirectory(src)
add_subdirectory(test)

include(CheckLanguage)
check_language(CUDA)

message("CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
message("CMAKE_CUDA_COMPILER_ID = ${CMAKE_CUDA_COMPILER_ID}")
message("CMAKE_CUDA_COMPILER_VERSION = ${CMAKE_CUDA_COMPILER_VERSION}")