cmake_minimum_required(VERSION 3.20)

message("configure gtest")

enable_testing()

find_package(GTest REQUIRED)

include(GoogleTest)

add_executable(
    test_fmm 
    test_fmm.cpp
)

target_include_directories(
    test_fmm PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(
    test_fmm
    core
    fftw3
    GTest::GTest
    #GTest::Main
)

gtest_discover_tests(test_fmm)
