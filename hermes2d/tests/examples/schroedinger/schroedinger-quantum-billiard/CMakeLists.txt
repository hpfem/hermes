IF(H2D_COMPLEX)

project(example-quantum-billiard)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-quantum-billiard ${BIN})
set_tests_properties(example-quantum-billiard PROPERTIES LABELS slow)

ENDIF(H2D_COMPLEX)
