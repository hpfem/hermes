project(test-calc-tetra-h1-neumann)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

# FIXME: THESE TESTS NEED TO BE FIXED AND ACTIVATED.
#add_test(${PROJECT_NAME}-1 ${BIN} tetra1-ref.mesh3d 2)
#add_test(${PROJECT_NAME}-2 ${BIN} tetra1.mesh3d 2)
