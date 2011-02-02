project(h3d-quadrature)

if(H3D_REAL)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_REAL})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(quadrature-1 ${BIN})

endif(H3D_REAL)
