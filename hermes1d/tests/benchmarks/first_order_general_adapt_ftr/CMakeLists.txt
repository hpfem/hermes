project(benchmark-first_order_general_adapt_ftr)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-first_order_general_adapt_ftr ${BIN})
set_tests_properties(benchmark-first_order_general_adapt_ftr PROPERTIES COST -500)
