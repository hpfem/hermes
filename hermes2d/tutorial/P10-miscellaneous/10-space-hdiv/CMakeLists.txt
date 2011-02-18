if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)
project(space-hdiv)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})