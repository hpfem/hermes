project(calc-hex-h1-mf)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

set(BINARIES) 

macro(ADD_BIN BIN DEFINE)
	add_executable(${BIN}
		hex-matfree.cpp
		${HERMES_COMMON_DIR}/timer.cpp
		${HERMES_COMMON_DIR}/error.cpp
		${HERMES_COMMON_DIR}/utils.cpp
	)
	
	set_target_properties(${BIN} PROPERTIES COMPILE_FLAGS -D${DEFINE})
	
	target_link_libraries(${BIN} ${TESTS_LIBRARIES})
	target_link_libraries(${BIN} ${HERMES_REAL_BIN})
	
	set(BINARIES ${BINARIES} ${BIN})
endmacro(ADD_BIN)

ADD_BIN(${PROJECT_NAME}-dir LIN_DIRICHLET)
ADD_BIN(${PROJECT_NAME}-neu LIN_NEUMANN)
ADD_BIN(${PROJECT_NAME}-nwt LIN_NEWTON)
ADD_BIN(${PROJECT_NAME}-nln NLN_DIRICHLET)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

foreach(B ${BINARIES})
	set(BIN ${PROJECT_BINARY_DIR}/${B})
	
#	add_test(${B}-1  ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d)
	add_test(${B}-4  ${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d)
#	add_test(${B}-27 ${BIN} ${MESHES_DIR}/mesh3d/hex27.mesh3d)
endforeach(B)

