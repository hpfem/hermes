project(tutorial-P04-01-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-P04-01-adapt ${BIN})
set_tests_properties(tutorial-P04-01-adapt PROPERTIES LABELS slow)
