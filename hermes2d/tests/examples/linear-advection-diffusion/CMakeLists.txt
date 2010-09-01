if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-linear-advection-diffusion)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-linear-advection-diffusion ${BIN})
set_tests_properties(example-linear-advection-diffusion PROPERTIES LABELS slow)
