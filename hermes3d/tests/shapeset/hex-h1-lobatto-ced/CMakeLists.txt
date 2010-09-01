project(hex-h1-lobatto-shapeset-ced)

iinclude(CMake.vars OPTIONAL)

if(WITH_HEX)
#
add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/timer.cpp
	${HERMES_COMMON_DIR}/error.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(hex-h1-lobbatto-shapeset-ced-1 ${BIN})
set_tests_properties(hex-h1-lobbatto-shapeset-ced-1 PROPERTIES LABELS long)

endif(WITH_HEX)

