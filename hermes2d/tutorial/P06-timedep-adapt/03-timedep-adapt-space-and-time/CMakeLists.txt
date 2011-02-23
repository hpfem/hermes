if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)

project(P06-03-timedep-adapt-space-and-time)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

add_subdirectory(tests)