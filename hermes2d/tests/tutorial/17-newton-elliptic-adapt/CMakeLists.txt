project(tutorial-17-newton-elliptic-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-17-newton-elliptic-adapt ${BIN})
set_tests_properties(tutorial-17-newton-elliptic-adapt PROPERTIES LABELS slow)
