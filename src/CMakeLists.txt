file(GLOB_RECURSE C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
list(REMOVE_ITEM C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.c")

set(PARSER_DIR ${PROJECT_SOURCE_DIR}/vendor/parser/src)

# for parser
file(GLOB PARSER_SOURCES
    ${PARSER_DIR}/parser.c
    ${PARSER_DIR}/symbol_list.c
    ${PARSER_DIR}/grammar.c
    ${PARSER_DIR}/memory.c
    ${PARSER_DIR}/result_manager.c
    ${PARSER_DIR}/rsm.c
    ${PARSER_DIR}/adapters/*.c
)

if (BUILD_SHARED_LIBS)
    add_library(cfl_lib SHARED
    ${C_SOURCES}
    ${PARSER_SOURCES}
)
else ()
    add_library(cfl_lib STATIC
    ${C_SOURCES}
    ${PARSER_SOURCES}
)
endif()

target_compile_options(cfl_lib PRIVATE -Wall -Wextra -Wpedantic)

target_include_directories(cfl_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PARSER_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/vendor/uthash>
    $<INSTALL_INTERFACE:include/cfl-idr-with-la>

)

set_target_properties(cfl_lib PROPERTIES OUTPUT_NAME "cfl-idr-with-la"
            PUBLIC_HEADER "include/cfl_idr.h" )
