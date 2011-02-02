project(lobatto-2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(lobatto-2 ${BIN})
set_tests_properties(lobatto-2 PROPERTIES COST -13.3)