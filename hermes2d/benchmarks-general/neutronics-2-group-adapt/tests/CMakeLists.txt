if(NOT H2D_REAL)
	return()
endif(NOT H2D_REAL)
project(test-benchmark-neutronics-2-group-adapt)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-benchmark-neutronics-2-group-adapt ${BIN})
set_tests_properties(test-benchmark-neutronics-2-group-adapt PROPERTIES LABELS slow)


