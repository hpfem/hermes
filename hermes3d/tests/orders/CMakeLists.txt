project(orders)

if(H3D_REAL)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/trace.cpp
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
add_test(orders-1 ${BIN})

endif(H3D_REAL)
