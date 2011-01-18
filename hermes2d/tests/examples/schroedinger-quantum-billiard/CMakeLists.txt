IF(H2D_COMPLEX)

project(example-quantum-billiard)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-quantum-billiard ${BIN})
set_tests_properties(example-quantum-billiard PROPERTIES LABELS slow)

ENDIF(H2D_COMPLEX)
