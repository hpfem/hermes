project(test-example-system_chaotic)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-system_chaotic ${BIN})
set_tests_properties(test-example-system_chaotic PROPERTIES COST -80)
