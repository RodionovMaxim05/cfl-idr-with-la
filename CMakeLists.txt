cmake_minimum_required(VERSION 3.20)

project(cfl-idr-with-la LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_EXECUTABLE "Build main executable" OFF)

set(LAGRAPH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LAGRAPH_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(LAGRAPH_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(src)

target_link_libraries(cfl_lib PUBLIC
    LAGraph
    LAGraphX
    GraphBLAS
)

if(BUILD_EXECUTABLE)
    add_executable(cfl-idr-with-la src/main.c)
    target_link_libraries(cfl-idr-with-la PRIVATE cfl_lib)
endif()
install(TARGETS cfl_lib
    EXPORT cfl-idr-with-laTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/cfl-idr-with-la
)

# Install public header
install(FILES include/cfl_idr.h DESTINATION include/cfl-idr-with-la)
# Export targets for installation
install(EXPORT cfl-idr-with-laTargets
    FILE cfl-idr-with-laTargets.cmake
    NAMESPACE cfl-idr-with-la::
    DESTINATION lib/cmake/cfl-idr-with-la
)
# Export targets in build tree (for FetchContent users)
export(TARGETS cfl_lib
    NAMESPACE cfl-idr-with-la::
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cfl-idr-with-laTargets.cmake
)
# Optional: Add alias for namespaced target
add_library(cfl-idr-with-la::cfl_lib ALIAS cfl_lib)

add_executable(wow src/main.c)
target_link_libraries(wow PRIVATE cfl_lib)
