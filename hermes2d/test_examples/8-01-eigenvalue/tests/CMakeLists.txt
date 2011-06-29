project(test-tutorial-P05-01-eigenvalue)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P05-01-eigenvalue ${BIN})
set_tests_properties(test-tutorial-P05-01-eigenvalue PROPERTIES LABELS slow)
