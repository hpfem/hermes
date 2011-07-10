project(test-cand_proj)

add_executable(${PROJECT_NAME} main.cpp functions.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-candidate_projection-quads ${BIN})
set_tests_properties(test-candidate_projection-quads PROPERTIES LABELS slow)
add_test(test-candidate_projection-tris ${BIN} -tri)
