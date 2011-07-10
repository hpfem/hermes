project(test-integrals-domain-perimeter)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-integrals-domain-perimeter-1 "${BIN}" domain.mesh0 8.0)
add_test(test-integrals-domain-perimeter-2 "${BIN}" domain.mesh1 6.650281539872986)
add_test(test-integrals-domain-perimeter-3 "${BIN}" domain.mesh2 10.681417125355473)
