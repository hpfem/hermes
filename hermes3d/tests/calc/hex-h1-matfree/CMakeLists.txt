project(calc-hex-h1-mf)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set(BINARIES) 

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}		hex-matfree.cpp)
	
  set_common_target_properties(${BIN})

	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS -D${DEFINE})
	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-dir LIN_DIRICHLET)
ADD_BIN(${PROJECT_NAME}-neu LIN_NEUMANN)
ADD_BIN(${PROJECT_NAME}-nwt LIN_NEWTON)
ADD_BIN(${PROJECT_NAME}-nln NLN_DIRICHLET)

# Tests

foreach(B ${BINARIES})
	set(BIN ${PROJECT_BINARY_DIR}/${B})
	
#	add_test(${B}-1  ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d)
	add_test(${B}-4  ${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d)
#	add_test(${B}-27 ${BIN} ${MESHES_DIR}/mesh3d/hex27.mesh3d)
endforeach(B)

