IF(H2D_COMPLEX)

project(test-example-waveguide)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-waveguide ${BIN})
set_tests_properties(test-example-waveguide PROPERTIES LABELS slow)

ENDIF(H2D_COMPLEX)
