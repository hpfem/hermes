project(test-adapt-exact-sin-H1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-adapt-exact-sin-H1 ${BIN})
set_tests_properties(test-adapt-exact-sin-H1 PROPERTIES COST -21.91)