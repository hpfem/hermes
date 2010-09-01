if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(tutorial-08-system)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-08-system ${BIN})

