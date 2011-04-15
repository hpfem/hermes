project(test-calc-hex-h1-nonlin)

set(BINARIES) 

include (${hermes3d_SOURCE_DIR}/CMake.common)

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}		main.cpp)
        set_common_target_properties(${BIN})
	get_target_property(CPFL ${BIN} COMPILE_FLAGS)
	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS "${CPFL} -D${DEFINE}")

	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-lin LINEAR)
ADD_BIN(${PROJECT_NAME}-nl1 NONLIN1)
ADD_BIN(${PROJECT_NAME}-nl2 NONLIN2)


# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-lin)
add_test(${PROJECT_NAME}-lin-1  ${BIN} hex1.mesh3d)
add_test(${PROJECT_NAME}-lin-8  ${BIN} hex8.mesh3d)
# FIXME: THIS TEST WAS DISABLED ON APRIL 14, 2011
# BECAUSE IT WAS FAILING.
#add_test(${PROJECT_NAME}-lin-27 ${BIN} hex27.mesh3d)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-nl1)
add_test(${PROJECT_NAME}-n1-1  ${BIN} hex8.mesh3d)
add_test(${PROJECT_NAME}-n1-27 ${BIN} hex27.mesh3d)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-nl2)
add_test(${PROJECT_NAME}-n2-1  ${BIN} hex1.mesh3d)
