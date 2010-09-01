if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)


project(complex-adapt-long)

set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
