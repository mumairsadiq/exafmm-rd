cmake_minimum_required(VERSION 3.20)

message("configure src")

find_package(BLAS REQUIRED)
message(STATUS BLAS found: ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
message(STATUS LAPACK found: ${LAPACK_LIBRARIES})

find_package(OpenMP REQUIRED)

add_library(core
    body.cpp
    surface.cpp
    fmm.cpp
    kernel.cpp
    mathfunc.cpp
    tree.cpp
    traverser.cpp
    ewald.cpp
    argument.cpp
    epkifmm.cpp
)

target_include_directories(
    core PRIVATE 
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(core PUBLIC
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    BLAS::BLAS
)

