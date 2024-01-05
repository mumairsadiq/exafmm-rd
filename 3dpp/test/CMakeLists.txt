cmake_minimum_required(VERSION 3.20)

message("configure test")

set(exec
    test_body
    test_surface
    test_blas
    test_lapack
    test_svd
    test_cell
    test_p2p
    test_p2m
    test_p2m2l2p
    test_fmm
    test_pbc
    test_p2m2m2l2p
    test_fftw
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
        fftw3
    )
endforeach()
