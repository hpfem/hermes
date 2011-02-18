project(adapt-exact-quadr-H1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-quadr-H1 ${BIN})
set_tests_properties(adapt-exact-quadr-H1 PROPERTIES COST -4.6)

