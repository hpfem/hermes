if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(test-example-euler-gamm-explicit)

add_executable(${PROJECT_NAME} main.cpp ../../euler_util.cpp ../../numerical_flux.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME}) 
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-euler-gamm-explicit ${BIN})
set_tests_properties(test-example-euler-gamm-explicit PROPERTIES LABELS slow)

