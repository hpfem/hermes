if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-optimal-meshes)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-optimal-meshes ${BIN} 0 1 2 3 4 10 10 10 10)
