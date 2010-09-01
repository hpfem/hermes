project(mesh-copy)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	mesh-copy.cpp
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

# copy
add_test(${PROJECT_NAME}-1-1 sh -c "${BIN} copy ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1-1.dump")
add_test(${PROJECT_NAME}-1-2 sh -c "${BIN} copy ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1-2.dump")
add_test(${PROJECT_NAME}-1-3 sh -c "${BIN} copy ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1-3.dump")

# copy base
add_test(${PROJECT_NAME}-base-1-1 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-1.dump")
add_test(${PROJECT_NAME}-base-1-2 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-1.dump")
add_test(${PROJECT_NAME}-base-1-3 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-1.dump")
add_test(${PROJECT_NAME}-base-1-4 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex1.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-1.dump")

add_test(${PROJECT_NAME}-base-2-1 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-2.dump")
add_test(${PROJECT_NAME}-base-2-2 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex2.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-2.dump")
add_test(${PROJECT_NAME}-base-2-3 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex2.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-2.dump")
add_test(${PROJECT_NAME}-base-2-4 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex2.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-2.dump")

add_test(${PROJECT_NAME}-base-4-1 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-4.dump")
add_test(${PROJECT_NAME}-base-4-2 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex4.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-4.dump")
add_test(${PROJECT_NAME}-base-4-3 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex4.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-4.dump")
add_test(${PROJECT_NAME}-base-4-4 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex4.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-4.dump")

add_test(${PROJECT_NAME}-base-8-1 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-8.dump")
add_test(${PROJECT_NAME}-base-8-2 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex8.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-8.dump")
add_test(${PROJECT_NAME}-base-8-3 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex8.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-8.dump")
add_test(${PROJECT_NAME}-base-8-4 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex8.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-8.dump")

add_test(${PROJECT_NAME}-base-27-1 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-27.dump")
add_test(${PROJECT_NAME}-base-27-2 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex27.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-27.dump")
add_test(${PROJECT_NAME}-base-27-3 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex27.mesh3d 1 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-27.dump")
add_test(${PROJECT_NAME}-base-27-4 sh -c "${BIN} copy_base ${MESHES_DIR}/mesh3d/hex27.mesh3d 1 z | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/base-27.dump")
