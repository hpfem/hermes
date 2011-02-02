project(multi-mesh)

include (${hermes3d_SOURCE_DIR}/CMake.common)

set(BINARIES)

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}		main.cpp)
	set_common_target_properties(${BIN})

	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS -D${DEFINE})

	target_link_libraries(${BIN} ${TESTS_LIBRARIES})
	target_link_libraries(${BIN} ${HERMES3D_REAL})

	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-rhs2 RHS2)
ADD_BIN(${PROJECT_NAME}-sys3 SYS3)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-rhs2)
add_test(${PROJECT_NAME}-rhs2-1  ${BIN} hex8.mesh3d)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-sys3)
add_test(${PROJECT_NAME}-sys3-1  ${BIN} hex8.mesh3d)
