project(tutorial-10-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-10-adapt ${BIN})
set_tests_properties(tutorial-10-adapt PROPERTIES LABELS slow)
