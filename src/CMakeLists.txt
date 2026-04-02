file(GLOB_RECURSE C_SOURCES ./*.c)
list(REMOVE_ITEM C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.c")

set(PARSER_DIR ${CMAKE_SOURCE_DIR}/vendor/parser/src)

add_library(cfl_lib STATIC
    ${C_SOURCES}
    ${PARSER_DIR}/parser.c
    ${PARSER_DIR}/symbol_list.c
)

set_target_properties(cfl_lib PROPERTIES OUTPUT_NAME "cfl-idr-with-la")

target_compile_options(cfl_lib PRIVATE -Wall -Wextra -Wpedantic) # -Werror)

target_include_directories(cfl_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PARSER_DIR}
    ${CMAKE_SOURCE_DIR}/vendor/lagraph/include/suitesparse
    ${CMAKE_SOURCE_DIR}/vendor/uthash
)
