project(example-first_order_general)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-first_order_general ${BIN})
set_tests_properties(example-first_order_general PROPERTIES COST -0.48)