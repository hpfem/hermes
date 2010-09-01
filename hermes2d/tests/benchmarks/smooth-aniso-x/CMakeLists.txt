if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-smooth-aniso-x)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-smooth-aniso-x-projection ${BIN})
add_test(benchmark-smooth-aniso-x-coarse ${BIN} -coarse_mesh)

