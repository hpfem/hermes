IF(H2D_COMPLEX)

project(bubbles)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(bubbles_quads ${BIN})

ENDIF(H2D_COMPLEX)
