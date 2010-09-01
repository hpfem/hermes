project(multi-mesh)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

set(BINARIES)

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}
		multi-mesh.cpp
		${HERMES_COMMON_DIR}/timer.cpp
		${HERMES_COMMON_DIR}/error.cpp
		${HERMES_COMMON_DIR}/utils.cpp
	)

	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS -D${DEFINE})

	target_link_libraries(${BIN} ${TESTS_LIBRARIES})
	target_link_libraries(${BIN} ${HERMES_REAL_BIN})

	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-rhs2 RHS2)
ADD_BIN(${PROJECT_NAME}-sys3 SYS3)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-rhs2)
add_test(${PROJECT_NAME}-rhs2-1  ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-sys3)
add_test(${PROJECT_NAME}-sys3-1  ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d)
