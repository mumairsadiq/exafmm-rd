cmake_minimum_required(VERSION 3.20)
project(rtfmm_3dp VERSION 1.0 LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 20)

message("configure all")

option(DEFINE_TEST_TREE "test tree" OFF)
if(DEFINE_TEST_TREE)
    add_compile_definitions(TEST_TREE)
endif()
message("DEFINE_TEST_TREE = ${DEFINE_TEST_TREE}")

add_subdirectory(src)
add_subdirectory(test)

include(CheckLanguage)
check_language(CUDA)

message("CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
message("CMAKE_CUDA_COMPILER_ID = ${CMAKE_CUDA_COMPILER_ID}")
message("CMAKE_CUDA_COMPILER_VERSION = ${CMAKE_CUDA_COMPILER_VERSION}")
message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
get_directory_property( DirDefs COMPILE_DEFINITIONS)
message( "COMPILE_DEFINITIONS = ${DirDefs}" )