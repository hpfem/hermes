project(adapt-exact-quadr-L2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-quadr-L2 ${BIN})
set_tests_properties(adapt-exact-quadr-L2 PROPERTIES COST -4.6)