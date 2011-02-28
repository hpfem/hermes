project(test-calc-hex-h1-mf)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set(BINARIES) 

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}		main.cpp)
	
    set_common_target_properties(${BIN})
    get_target_property(CPFL ${BIN} COMPILE_FLAGS)
	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS "${CPFL} -D${DEFINE}")
	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-dir LIN_DIRICHLET)
ADD_BIN(${PROJECT_NAME}-neu LIN_NEUMANN)
ADD_BIN(${PROJECT_NAME}-nwt LIN_NEWTON)
ADD_BIN(${PROJECT_NAME}-nln NLN_DIRICHLET)

# Tests

foreach(B ${BINARIES})
	set(BIN ${PROJECT_BINARY_DIR}/${B})
	
#	add_test(${B}-1  ${BIN} hex1.mesh3d)
	add_test(${B}-4  ${BIN} hex4.mesh3d)
#	add_test(${B}-27 ${BIN} hex27.mesh3d)
endforeach(B)

