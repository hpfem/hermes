if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-lshape)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-lshape-projection ${BIN})
add_test(benchmark-lshape-coarse ${BIN} -coarse_mesh)
set_tests_properties(benchmark-lshape-projection PROPERTIES LABELS slow)
set_tests_properties(benchmark-lshape-coarse PROPERTIES LABELS slow)
