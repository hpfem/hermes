if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(19-rk-comparison)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})