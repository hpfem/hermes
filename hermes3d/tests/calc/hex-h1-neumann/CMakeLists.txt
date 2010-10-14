project(calc-hex-h1-neumann)
add_executable(${PROJECT_NAME}	hex-neumann.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1    ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-sc-1 ${BIN} ${MESHES_DIR}/mesh3d/hex1-sc.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-2x   ${BIN} ${MESHES_DIR}/mesh3d/hex2x.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-2y   ${BIN} ${MESHES_DIR}/mesh3d/hex2y.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-2z   ${BIN} ${MESHES_DIR}/mesh3d/hex2z.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-4x   ${BIN} ${MESHES_DIR}/mesh3d/hex4x.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-4y   ${BIN} ${MESHES_DIR}/mesh3d/hex4y.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-4z   ${BIN} ${MESHES_DIR}/mesh3d/hex4z.mesh3d 2 2 2)
add_test(${PROJECT_NAME}-8    ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d  2 2 2)
