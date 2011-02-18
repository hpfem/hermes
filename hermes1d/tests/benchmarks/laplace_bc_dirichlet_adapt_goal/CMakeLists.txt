project(benchmark-laplace_bc_dirichlet_adapt_goal)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-laplace_bc_dirichlet_adapt_goal ${BIN})
set_tests_properties(benchmark-laplace_bc_dirichlet_adapt_goal PROPERTIES COST -10)