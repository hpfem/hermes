project(test-mesh_refinements)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-mesh_refinements-1 "${BIN}" domain.mesh)
add_test(test-mesh_refinements-2 "${BIN}" bracket.mesh)
add_test(test-mesh_refinements-3 "${BIN}" square.mesh)
add_test(test-mesh_refinements-4 "${BIN}" square_tri.mesh)
