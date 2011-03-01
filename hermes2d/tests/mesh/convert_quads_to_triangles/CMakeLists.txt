project(test-convert_quads_to_triangles)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-convert_quads_to_triangles-1 "${BIN}" domain.mesh)
add_test(test-convert_quads_to_triangles-2 "${BIN}" bracket.mesh)
add_test(test-convert_quads_to_triangles-3 "${BIN}" square.mesh)
add_test(test-convert_quads_to_triangles-4 "${BIN}" square_tri.mesh)

