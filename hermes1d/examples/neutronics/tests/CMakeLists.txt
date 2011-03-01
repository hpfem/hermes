project(test-example-neutronics)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-neutronics ${BIN})
#set_tests_properties(test-example-neutronics PROPERTIES WILL_FAIL TRUE)