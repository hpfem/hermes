if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-nist-05)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-nist-05 ${BIN})
set_tests_properties(benchmark-nist-05 PROPERTIES LABELS slow)
