project(convert_triangles_to_quads)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(convert_triangles_to_quads-1 "${BIN}" square.mesh)
add_test(convert_triangles_to_quads-2 "${BIN}" square_tri.mesh)
add_test(convert_triangles_to_quads-3 "${BIN}" domain.mesh)
add_test(convert_triangles_to_quads-4 "${BIN}" bracket.mesh)

