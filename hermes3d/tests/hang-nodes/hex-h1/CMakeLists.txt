project(hnnd-hex-h1)
set(BIN_DIRICHLET ${PROJECT_NAME}-dirichlet)
set(BIN_NEWTON ${PROJECT_NAME}-newton)

include (${hermes3d_SOURCE_DIR}/CMake.common)

# Dirichlet ####

add_executable(${BIN_DIRICHLET}	main.cpp)
set_common_target_properties(${BIN_DIRICHLET})

set_target_properties(${BIN_DIRICHLET}
    PROPERTIES
	COMPILE_FLAGS -DDIRICHLET)

# Newton ####

add_executable(${BIN_NEWTON}	main.cpp)
set_common_target_properties(${BIN_NEWTON})

set_target_properties(${BIN_NEWTON}
    PROPERTIES
	COMPILE_FLAGS -DNEWTON)
	
# Tests

set(BIN ${PROJECT_BINARY_DIR}/${BIN_DIRICHLET})
set(NAME ${PROJECT_NAME}-dir)
include(../tests-hex.cmake)

set(BIN ${PROJECT_BINARY_DIR}/${BIN_NEWTON})
set(NAME ${PROJECT_NAME}-nwt)
include(../tests-hex.cmake)
