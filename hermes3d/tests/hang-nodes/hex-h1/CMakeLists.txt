project(hnnd-hex-h1)

include(CMake.vars OPTIONAL)

set(BIN_DIRICHLET ${PROJECT_NAME}-dirichlet)
set(BIN_NEWTON ${PROJECT_NAME}-newton)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

# Dirichlet ####

add_executable(${BIN_DIRICHLET}
	hang-nodes-cont.cpp
	${HERMES_COMMON_DIR}/timer.cpp
	${HERMES_COMMON_DIR}/error.cpp
)

set_target_properties(${BIN_DIRICHLET}
    PROPERTIES
	COMPILE_FLAGS -DDIRICHLET)

target_link_libraries(${BIN_DIRICHLET} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_DIRICHLET} ${HERMES_REAL_BIN})

# Newton ####

add_executable(${BIN_NEWTON}
	hang-nodes-cont.cpp
	${HERMES_COMMON_DIR}/timer.cpp
	${HERMES_COMMON_DIR}/error.cpp
)

set_target_properties(${BIN_NEWTON}
    PROPERTIES
	COMPILE_FLAGS -DNEWTON)

target_link_libraries(${BIN_NEWTON} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_NEWTON} ${HERMES_REAL_BIN})


configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${BIN_DIRICHLET})
set(NAME ${PROJECT_NAME}-dir)
include(../tests-hex.cmake)

set(BIN ${PROJECT_BINARY_DIR}/${BIN_NEWTON})
set(NAME ${PROJECT_NAME}-nwt)
include(../tests-hex.cmake)
