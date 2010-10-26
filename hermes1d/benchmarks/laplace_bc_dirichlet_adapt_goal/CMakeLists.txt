project(laplace_bc_dirichlet_adapt_goal)

add_executable(${PROJECT_NAME} main.cpp)
include(../CMake.common)

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})