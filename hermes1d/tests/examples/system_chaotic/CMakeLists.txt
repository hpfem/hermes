project(example-system_chaotic)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-system_chaotic ${BIN})
set_tests_properties(example-system_chaotic PROPERTIES COST -80)
set_tests_properties(example-system_chaotic PROPERTIES WILL_FAIL TRUE)