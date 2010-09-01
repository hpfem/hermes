IF(H2D_COMPLEX)

project(bubbles)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(bubbles_quads ${BIN})

ENDIF(H2D_COMPLEX)
