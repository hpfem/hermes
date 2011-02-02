project(example-laplace_bc_newton2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-laplace_bc_newton2 ${BIN})
set_tests_properties(example-laplace_bc_newton2 PROPERTIES WILL_FAIL TRUE)