project(example-system_exp)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-system_exp ${BIN})
set_tests_properties(example-system_exp PROPERTIES WILL_FAIL TRUE)