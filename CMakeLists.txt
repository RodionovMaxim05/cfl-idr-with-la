cmake_minimum_required(VERSION 3.20)

project(cfl-idr-with-la LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(LAGRAPH_DIR "" CACHE PATH "Path to LAGraph source directory")
if(NOT LAGRAPH_DIR)
    message(FATAL_ERROR "Please set LAGRAPH_DIR: cmake -DLAGRAPH_DIR=/path/to/lagraph ..")
endif()

set(LAGRAPH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LAGRAPH_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(LAGRAPH_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(${LAGRAPH_DIR} lagraph_build)

add_subdirectory(src)

target_link_libraries(cfl_lib PRIVATE
    LAGraph
    LAGraphX
    /usr/local/lib/libgraphblas.so
)

add_executable(cfl-idr-with-la src/main.c)
target_link_libraries(cfl-idr-with-la PRIVATE cfl_lib)

include(CTest)
add_subdirectory(tests)
