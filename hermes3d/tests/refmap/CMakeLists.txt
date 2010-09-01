project(ref-map)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/trace.cpp
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

add_test(${PROJECT_NAME}-1 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex1.dump")
add_test(${PROJECT_NAME}-2 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex2.dump")
add_test(${PROJECT_NAME}-4 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex4.dump")
add_test(${PROJECT_NAME}-8 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex8.dump")
add_test(${PROJECT_NAME}-27 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex27.dump")

