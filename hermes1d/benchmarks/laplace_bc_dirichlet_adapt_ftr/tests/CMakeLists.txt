project(test-benchmark-laplace_bc_dirichlet_adapt_ftr)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-benchmark-laplace_bc_dirichlet_adapt_ftr ${BIN})
set_tests_properties(test-benchmark-laplace_bc_dirichlet_adapt_ftr PROPERTIES COST -10)