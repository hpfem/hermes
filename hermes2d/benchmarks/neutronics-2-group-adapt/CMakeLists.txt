if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)
project(neutronics-2-group-adapt)

add_executable(${PROJECT_NAME} main.cpp scaled_scalar_view.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})