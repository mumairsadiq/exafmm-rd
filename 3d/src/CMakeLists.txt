cmake_minimum_required(VERSION 3.20)

message("configure src")

add_library(core
    body.cpp
    surface.cpp
)

target_include_directories(
    core PRIVATE 
    "${PROJECT_SOURCE_DIR}/include"
)
