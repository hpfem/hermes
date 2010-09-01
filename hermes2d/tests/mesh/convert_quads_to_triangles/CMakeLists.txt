project(convert_quads_to_triangles)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(convert_quads_to_triangles-1 "${BIN}" domain.mesh)
add_test(convert_quads_to_triangles-2 "${BIN}" bracket.mesh)
add_test(convert_quads_to_triangles-3 "${BIN}" square.mesh)
add_test(convert_quads_to_triangles-4 "${BIN}" square_tri.mesh)

