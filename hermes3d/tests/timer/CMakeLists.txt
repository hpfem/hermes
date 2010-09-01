project(timer)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${hermes3d_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/trace.cpp
	${HERMES_COMMON_DIR}/timer.cpp
)

target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

if(H3D_REAL)
	target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
elseif(H3D_COMPLEX)
	target_link_libraries(${PROJECT_NAME} ${HERMES_CPLX_BIN})
endif(H3D_REAL)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(timer-1 ${BIN})

