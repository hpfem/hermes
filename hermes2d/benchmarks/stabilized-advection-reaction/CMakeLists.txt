if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(stabilized-advection-reaction)

add_executable(${PROJECT_NAME} main.cpp exact.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})