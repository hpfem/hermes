if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-np)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-np ${BIN})
set_tests_properties(example-np PROPERTIES LABELS slow)
