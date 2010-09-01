project(shapeset-hex-hcurl)

include(CMake.vars OPTIONAL)

if(WITH_HEX)
#
add_executable(${PROJECT_NAME}
	main.cpp
	lindep.cpp
	zero.cpp
	cont.cpp
	grad.cpp
	graddir.cpp
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
add_test(${PROJECT_NAME}-1 ${BIN})
set_tests_properties(${PROJECT_NAME}-1 PROPERTIES LABELS long)

endif(WITH_HEX)
