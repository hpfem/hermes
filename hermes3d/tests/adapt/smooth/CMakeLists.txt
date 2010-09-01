project(adapt-smooth)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/timer.cpp
	${HERMES_COMMON_DIR}/error.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

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
