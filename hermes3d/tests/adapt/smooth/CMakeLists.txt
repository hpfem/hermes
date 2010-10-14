project(adapt-smooth)
add_executable(${PROJECT_NAME}	main.cpp)

include(${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-x ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d x)
add_test(${PROJECT_NAME}-y ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d y)
add_test(${PROJECT_NAME}-z ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d z)
add_test(${PROJECT_NAME}-xy ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d xy)
add_test(${PROJECT_NAME}-xz ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d xz)
add_test(${PROJECT_NAME}-yz ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d yz)
add_test(${PROJECT_NAME}-xyz ${BIN} ${MESHES_DIR}/mesh3d/hex-0-1.mesh3d xyz)

set_tests_properties(${PROJECT_NAME}-x PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-y PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-z PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-xy PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-xz PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-yz PROPERTIES LABELS long)
set_tests_properties(${PROJECT_NAME}-xyz PROPERTIES LABELS long)
