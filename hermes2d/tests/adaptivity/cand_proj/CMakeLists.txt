project(cand_proj)

add_executable(${PROJECT_NAME} main.cpp functions.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(candidate_projection-quads ${BIN})
set_tests_properties(candidate_projection-quads PROPERTIES LABELS slow)
add_test(candidate_projection-tris ${BIN} -tri)
