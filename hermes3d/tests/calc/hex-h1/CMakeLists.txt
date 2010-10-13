project(calc-hex-h1)
add_executable(${PROJECT_NAME}	hex.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-x90 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-x90.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-y90 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-y90.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-z90 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-z90.mesh3d 2)
add_test(${PROJECT_NAME}-2  ${BIN} ${MESHES_DIR}/mesh3d/hex2.mesh3d 2)
add_test(${PROJECT_NAME}-2x ${BIN} ${MESHES_DIR}/mesh3d/hex2x.mesh3d 2)
add_test(${PROJECT_NAME}-2y ${BIN} ${MESHES_DIR}/mesh3d/hex2y.mesh3d 2)
add_test(${PROJECT_NAME}-2z ${BIN} ${MESHES_DIR}/mesh3d/hex2z.mesh3d 2)
add_test(${PROJECT_NAME}-4  ${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d 2)
add_test(${PROJECT_NAME}-4x ${BIN} ${MESHES_DIR}/mesh3d/hex4x.mesh3d 2)
add_test(${PROJECT_NAME}-4y ${BIN} ${MESHES_DIR}/mesh3d/hex4y.mesh3d 2)
add_test(${PROJECT_NAME}-4z ${BIN} ${MESHES_DIR}/mesh3d/hex4z.mesh3d 2)
add_test(${PROJECT_NAME}-8  ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d 2)
