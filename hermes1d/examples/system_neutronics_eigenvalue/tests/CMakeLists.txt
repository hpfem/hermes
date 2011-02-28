project(test-example-system_neutronics_eigenvalue)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-system_neutronics_eigenvalue ${BIN})
#set_tests_properties(test-example-system_neutronics_eigenvalue PROPERTIES WILL_FAIL TRUE)