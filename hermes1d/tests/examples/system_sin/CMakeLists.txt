project(example-system_sin)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-system_sin ${BIN})
set_tests_properties(example-system_sin PROPERTIES COST -0.00224)
