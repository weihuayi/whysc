cmake_minimum_required(VERSION 3.0)

project(WHYSC)

include_directories(${PROJECT_SOURCE_DIR}/include)
message(STATUS "${PROJECT_SOURCE_DIR}")

add_subdirectory(test)
