project(hex-h1-1face)

if(DEV_TESTS)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	main.cpp
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

enable_testing()
add_test(${PROJECT_NAME}-1 ${PROJECT_NAME})

endif(DEV_TESTS)
