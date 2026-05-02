file(GLOB_RECURSE C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
list(REMOVE_ITEM C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.c")

set(PARSER_DIR ${CMAKE_SOURCE_DIR}/vendor/parser/src)

# for parser
file(GLOB PARSER_SOURCES
    ${PARSER_DIR}/parser.c
    ${PARSER_DIR}/symbol_list.c
    ${PARSER_DIR}/grammar.c
    ${PARSER_DIR}/memory.c
    ${PARSER_DIR}/result_manager.c
    ${PARSER_DIR}/adapters/*.c
)

add_library(cfl_lib STATIC
    ${C_SOURCES}
    ${PARSER_SOURCES}
)

set_target_properties(cfl_lib PROPERTIES OUTPUT_NAME "cfl-idr-with-la")

target_compile_options(cfl_lib PRIVATE -Wall -Wextra -Wpedantic) # -Werror)

target_include_directories(cfl_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PARSER_DIR}
    ${LAGRAPH_DIR}/include
    ${CMAKE_SOURCE_DIR}/vendor/uthash
    /usr/local/include/suitesparse
)
