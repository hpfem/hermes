project(judy-templates)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/trace.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${JUDY_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY}) 

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(judy-templates-1 ${BIN})

