project(test-example-system_pendulum)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-system_pendulum ${BIN})
set_tests_properties(test-example-system_pendulum PROPERTIES COST -22.10)