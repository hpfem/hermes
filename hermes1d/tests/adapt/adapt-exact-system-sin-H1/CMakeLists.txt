project(adapt-exact-system-sin-H1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-system-sin-H1 ${BIN})
set_tests_properties(adapt-exact-system-sin-H1 PROPERTIES COST -5.6)