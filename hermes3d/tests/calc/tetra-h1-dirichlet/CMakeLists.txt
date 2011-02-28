project(test-calc-tetra-h1-dirichlet)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})


# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} tetra1.mesh3d 4)
add_test(${PROJECT_NAME}-2 ${BIN} tetra2.mesh3d 3)
add_test(${PROJECT_NAME}-3 ${BIN} tetra4.mesh3d 2)
add_test(${PROJECT_NAME}-4 ${BIN} tetra8.mesh3d 2)
