cmake_minimum_required(VERSION 3.20)
project(rtfmm_3dirpp VERSION 1.0 LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 20)

#add_compile_options(-march=native)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2") 
#add_compile_options(-std=c++11 -fopenmp -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -msse4a -msha -maes -mavx -mfma -mavx2 -mrdrnd -mbmi -mbmi2 -madx -mabm -g -O2)

message("configure all")

option(DEFINE_TEST_TREE "test tree" OFF)
if(DEFINE_TEST_TREE)
    add_compile_definitions(TEST_TREE)
endif()
message("DEFINE_TEST_TREE = ${DEFINE_TEST_TREE}")

add_subdirectory(src)
add_subdirectory(test)


message ("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
get_directory_property( DirDefs COMPILE_DEFINITIONS)
message( "COMPILE_DEFINITIONS = ${DirDefs}" )