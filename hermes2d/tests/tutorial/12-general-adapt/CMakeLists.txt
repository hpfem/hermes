project(tutorial-12-general-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-12-general-adapt ${BIN})
set_tests_properties(tutorial-12-general-adapt PROPERTIES LABELS slow)

