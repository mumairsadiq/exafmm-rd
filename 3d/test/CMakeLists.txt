cmake_minimum_required(VERSION 3.20)

message("configure test")

set(exec
    test_body
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
    target_link_libraries(${filename} core)
endforeach()
