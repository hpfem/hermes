if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-singular-perturbation)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-singular-perturbation ${BIN})
set_tests_properties(example-singular-perturbation PROPERTIES LABELS slow)
