if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-sdirk-22)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-sdirk-22-e ${BIN} e)
add_test(benchmark-sdirk-22-s ${BIN} s)
