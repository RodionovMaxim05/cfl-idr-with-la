cmake_minimum_required(VERSION 3.20)

project(cfl-idr-with-la LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(LAGRAPH_DIR ${CMAKE_SOURCE_DIR}/vendor/lagraph)

add_subdirectory(src)

target_include_directories(cfl_lib PUBLIC
    ${LAGRAPH_DIR}/include/suitesparse
    /usr/local/include/suitesparse
)

target_link_libraries(cfl_lib PRIVATE
    ${LAGRAPH_DIR}/lib/liblagraph.so.1.2.1
    ${LAGRAPH_DIR}/lib/liblagraphx.so.1.2.1
    /usr/local/lib/libgraphblas.so
)

add_executable(cfl-idr-with-la src/main.c)
target_link_libraries(cfl-idr-with-la PRIVATE cfl_lib)

enable_testing()
add_subdirectory(tests)
