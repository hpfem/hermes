project(test-cand_proj)

add_executable(${PROJECT_NAME} main.cpp functions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-candidate_projection-quads ${BIN})
set_tests_properties(test-candidate_projection-quads PROPERTIES LABELS slow)
add_test(test-candidate_projection-tris ${BIN} -tri)
