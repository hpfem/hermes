IF(H3D_COMPLEX)

project(test-example-rough)
add_executable(${PROJECT_NAME} main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_CPLX})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME} ${BIN})

ENDIF(H3D_COMPLEX)
