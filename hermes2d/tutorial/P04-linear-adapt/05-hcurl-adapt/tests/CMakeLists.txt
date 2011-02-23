if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)
project(tutorial-P04-05-hcurl-adapt)

add_executable(${PROJECT_NAME} main.cpp bessel.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-P04-05-hcurl-adapt ${BIN})
set_tests_properties(tutorial-P04-05-hcurl-adapt PROPERTIES LABELS slow)

