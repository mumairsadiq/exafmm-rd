cmake_minimum_required(VERSION 3.20)

message("configure test")

find_package(BLAS REQUIRED)
message(STATUS BLAS found: ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
message(STATUS LAPACK found: ${LAPACK_LIBRARIES})

find_package(OpenMP REQUIRED)

set(exec
    test_body
    test_surface
    test_blas
    test_lapack
    test_svd
    test_cell
)

foreach(filename IN LISTS exec)
    message("configure ${filename}")
    add_executable(
        ${filename} 
        "${filename}.cpp"
    )
    target_include_directories(
        ${filename} PUBLIC 
        "${PROJECT_SOURCE_DIR}/include"
    )
    target_link_libraries(${filename} 
        core 
        OpenMP::OpenMP_CXX
        LAPACK::LAPACK
        BLAS::BLAS
    )
endforeach()
