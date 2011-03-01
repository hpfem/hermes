project(test-calc-hex-h1)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} hex1.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-x90 ${BIN} hex1-rot-x90.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-y90 ${BIN} hex1-rot-y90.mesh3d 2)
add_test(${PROJECT_NAME}-1-rot-z90 ${BIN} hex1-rot-z90.mesh3d 2)
add_test(${PROJECT_NAME}-2  ${BIN} hex2.mesh3d 2)
add_test(${PROJECT_NAME}-2x ${BIN} hex2x.mesh3d 2)
add_test(${PROJECT_NAME}-2y ${BIN} hex2y.mesh3d 2)
add_test(${PROJECT_NAME}-2z ${BIN} hex2z.mesh3d 2)
add_test(${PROJECT_NAME}-4  ${BIN} hex4.mesh3d 2)
add_test(${PROJECT_NAME}-4x ${BIN} hex4x.mesh3d 2)
add_test(${PROJECT_NAME}-4y ${BIN} hex4y.mesh3d 2)
add_test(${PROJECT_NAME}-4z ${BIN} hex4z.mesh3d 2)
add_test(${PROJECT_NAME}-8  ${BIN} hex8.mesh3d 2)
