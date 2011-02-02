if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)
project(tutorial-14-hcurl-adapt)

add_executable(${PROJECT_NAME} main.cpp bessel.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-14-hcurl-adapt ${BIN})
set_tests_properties(tutorial-14-hcurl-adapt PROPERTIES LABELS slow)

