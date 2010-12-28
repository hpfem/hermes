project(tutorial-15-newton-elliptic-1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-15-newton-elliptic-1 ${BIN})
