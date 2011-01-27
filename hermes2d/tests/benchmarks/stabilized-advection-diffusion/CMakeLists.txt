if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-stabilized-advection-diffusion)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-stabilized-advection-diffusion ${BIN})
set_tests_properties(benchmark-stabilized-advection-diffusion PROPERTIES LABELS slow)
