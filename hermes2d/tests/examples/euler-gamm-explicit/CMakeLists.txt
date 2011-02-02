if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-euler-gamm-explicit)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)
