project(example-homogeneous_line)
add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-homogeneous_line ${BIN})
set_tests_properties(example-homogeneous_line PROPERTIES COST -124.55)
set_tests_properties(example-homogeneous_line PROPERTIES WILL_FAIL TRUE)