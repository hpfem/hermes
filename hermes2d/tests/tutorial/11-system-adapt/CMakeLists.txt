project(tutorial-11-system-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-11-system-adapt ${BIN})
set_tests_properties(tutorial-11-system-adapt PROPERTIES LABELS slow)
