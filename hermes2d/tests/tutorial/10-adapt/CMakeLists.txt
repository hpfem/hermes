project(tutorial-10-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-10-adapt ${BIN})
set_tests_properties(tutorial-10-adapt PROPERTIES LABELS slow)
