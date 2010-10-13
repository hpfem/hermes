project(projection)
add_executable(${PROJECT_NAME}	main.cpp)

include(${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

#
# Tests
#

enable_testing()

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(projection-h1-1 ${BIN} h1 ${MESHES_DIR}/mesh3d/hex1.mesh3d)
add_test(projection-h1-2 ${BIN} h1 ${MESHES_DIR}/mesh3d/hex2.mesh3d)
add_test(projection-h1-3 ${BIN} h1 ${MESHES_DIR}/mesh3d/hex4.mesh3d)
add_test(projection-h1-4 ${BIN} h1 ${MESHES_DIR}/mesh3d/hex8.mesh3d)
add_test(projection-h1-5 ${BIN} h1 ${MESHES_DIR}/mesh3d/hex27.mesh3d)

add_test(projection-h1i-1 ${BIN} h1-ipol ${MESHES_DIR}/mesh3d/hex1.mesh3d)
add_test(projection-h1i-2 ${BIN} h1-ipol ${MESHES_DIR}/mesh3d/hex2.mesh3d)
add_test(projection-h1i-3 ${BIN} h1-ipol ${MESHES_DIR}/mesh3d/hex4.mesh3d)
add_test(projection-h1i-4 ${BIN} h1-ipol ${MESHES_DIR}/mesh3d/hex8.mesh3d)
add_test(projection-h1i-5 ${BIN} h1-ipol ${MESHES_DIR}/mesh3d/hex27.mesh3d)

set_tests_properties(projection-h1-5 PROPERTIES LABELS long)
set_tests_properties(projection-h1i-4 PROPERTIES LABELS long)
set_tests_properties(projection-h1i-5 PROPERTIES LABELS long)
