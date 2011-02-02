project(benchmark-system_sin_adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-system_sin_adapt ${BIN})
set_tests_properties(benchmark-system_sin_adapt PROPERTIES COST -9.99)