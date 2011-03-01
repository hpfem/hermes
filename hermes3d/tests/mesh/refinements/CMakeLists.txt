project(test-mesh-refinements)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1-1 "${BIN}" hex1.mesh3d 1 x)
add_test(${PROJECT_NAME}-1-2 "${BIN}" hex1.mesh3d 1 x 2 y)
add_test(${PROJECT_NAME}-4-1 "${BIN}" hex4.mesh3d 1 xyz)
