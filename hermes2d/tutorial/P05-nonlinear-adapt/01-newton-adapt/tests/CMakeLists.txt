project(tutorial-P05-01-newton-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-P05-01-newton-adapt ${BIN})
set_tests_properties(tutorial-P05-01-newton-adapt PROPERTIES LABELS slow)
