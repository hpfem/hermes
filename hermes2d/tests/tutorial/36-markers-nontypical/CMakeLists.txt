project(tutorial-36-markers-nontypical)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-36-markers-nontypical ${BIN})
