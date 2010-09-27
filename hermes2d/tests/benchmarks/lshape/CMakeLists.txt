if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-lshape)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-lshape ${BIN})
set_tests_properties(benchmark-lshape PROPERTIES LABELS slow)
