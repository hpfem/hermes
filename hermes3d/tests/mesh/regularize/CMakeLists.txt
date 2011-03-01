project(test-regularize)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1 "${BIN}" hex4.mesh3d 2 xyz)
add_test(${PROJECT_NAME}-2 "${BIN}" hex4.mesh3d 2 xyz 5 xyz)
