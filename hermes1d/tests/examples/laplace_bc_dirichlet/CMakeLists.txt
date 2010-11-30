project(example-laplace_bc_dirichlet)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-laplace_bc_dirichlet ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
