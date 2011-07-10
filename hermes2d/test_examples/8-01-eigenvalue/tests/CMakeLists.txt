project(test-tutorial-P05-01-eigenvalue)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
	
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P05-01-eigenvalue ${BIN})
set_tests_properties(test-tutorial-P05-01-eigenvalue PROPERTIES LABELS slow)
