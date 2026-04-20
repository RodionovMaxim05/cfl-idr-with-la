file(GLOB_RECURSE TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/test_*.c")

set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test_data/correctness")

foreach(TEST_SOURCE ${TEST_SOURCES})
    file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_SOURCE_DIR}" "${TEST_SOURCE}")

    string(REPLACE ".c" "" TEST_PATH "${REL_PATH}")

    string(REPLACE "/" "_" TEST_NAME "${TEST_PATH}")

    add_executable(${TEST_NAME} ${TEST_SOURCE})

    target_compile_definitions(${TEST_NAME} PRIVATE TEST_DATA_DIR="${TEST_DATA_DIR}")

    target_include_directories(${TEST_NAME} PRIVATE
        ${LAGRAPH_DIR}/include
        /usr/local/include/suitesparse
    )

    target_link_libraries(${TEST_NAME} PRIVATE cfl_lib)

    add_test(NAME ${TEST_NAME} COMMAND $<TARGET_FILE:${TEST_NAME}>)
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "cfl-idr")

    message(STATUS "Registered test: ${TEST_NAME} from ${REL_PATH}")
endforeach()
