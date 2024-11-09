cmake_minimum_required(VERSION 3.20)

message("configure test")

file(GLOB exec *.cpp)

foreach(file ${exec})
    get_filename_component(filename ${file} NAME_WE)
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
