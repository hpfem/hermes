project(test-hnnd-hex-h1)
set(BIN_DIRICHLET ${PROJECT_NAME}-dirichlet)
set(BIN_NEWTON ${PROJECT_NAME}-newton)

include (${hermes3d_SOURCE_DIR}/CMake.common)

# Dirichlet ####

add_executable(${BIN_DIRICHLET}	main.cpp)
set_common_target_properties(${BIN_DIRICHLET})
get_target_property(CPFL ${BIN_DIRICHLET} COMPILE_FLAGS)
set_target_properties(${BIN_DIRICHLET}
    PROPERTIES
	COMPILE_FLAGS "${CPFL} -DDIRICHLET")

# Newton ####

add_executable(${BIN_NEWTON}	main.cpp)
set_common_target_properties(${BIN_NEWTON})
get_target_property(CPFL ${BIN_NEWTON} COMPILE_FLAGS)
set_target_properties(${BIN_NEWTON}
    PROPERTIES
	COMPILE_FLAGS "${CPFL} -DNEWTON")
	
# Tests

set(BIN ${PROJECT_BINARY_DIR}/${BIN_DIRICHLET})
set(NAME ${PROJECT_NAME}-dir)
include(../tests-hex.cmake)

set(BIN ${PROJECT_BINARY_DIR}/${BIN_NEWTON})
set(NAME ${PROJECT_NAME}-nwt)
include(../tests-hex.cmake)
