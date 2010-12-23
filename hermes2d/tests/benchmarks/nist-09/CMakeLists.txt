if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-nist-9)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-nist-9 ${BIN})
set_tests_properties(benchmark-nist-9 PROPERTIES LABELS slow)
