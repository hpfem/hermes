if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(neutronics-heat-conduction-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
