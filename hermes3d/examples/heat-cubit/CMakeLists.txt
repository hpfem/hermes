if(WITH_EXODUSII)

project(heat-cubit)

include(CMake.vars OPTIONAL)


add_executable(${PROJECT_NAME}
	main.cpp
)

include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

endif(WITH_EXODUSII)
