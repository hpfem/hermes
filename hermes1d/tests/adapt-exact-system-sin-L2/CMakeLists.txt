project(adapt-exact-system-sin-L2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-system-sin-L2 ${BIN})
