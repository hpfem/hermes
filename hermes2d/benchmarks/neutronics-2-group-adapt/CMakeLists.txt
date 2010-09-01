if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)

project(neutronics-2-group-adapt)
set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp scaled_scalar_view.cpp)
include (../CMake.common)

