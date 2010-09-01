project(adapt-exact-quadr-H1-solver0)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-quadr-H1-solver0 ${BIN})
