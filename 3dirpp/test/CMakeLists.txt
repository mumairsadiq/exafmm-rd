cmake_minimum_required(VERSION 3.20)

message("configure test")

set(exec
    test_tree
    test_traverser
    test_fft
    test_kernel
    test_p2p
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
