project(hnnd-hex-hcurl)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	hang-nodes-cont.cpp
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
set(NAME ${PROJECT_NAME})
include(../tests-hex.cmake)

set_tests_properties(hnnd-hex-hcurl-1-01 PROPERTIES LABELS long)
set_tests_properties(hnnd-hex-hcurl-1-02 PROPERTIES LABELS long)
