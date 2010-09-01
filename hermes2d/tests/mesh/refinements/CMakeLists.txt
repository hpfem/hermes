project(refinements)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(refinements-1 "${BIN}" domain.mesh)
add_test(refinements-2 "${BIN}" bracket.mesh)
add_test(refinements-3 "${BIN}" square.mesh)
add_test(refinements-4 "${BIN}" square_tri.mesh)
