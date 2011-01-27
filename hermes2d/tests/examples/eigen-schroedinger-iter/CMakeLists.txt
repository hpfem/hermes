if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-eigen-schroedinger-iter)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-eigen-schroedinger-iter ${BIN})
