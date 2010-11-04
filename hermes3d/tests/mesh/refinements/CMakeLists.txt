project(mesh-refinements)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1-1 sh -c "${BIN} hex1.mesh3d 1 x | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1-1.dump")
add_test(${PROJECT_NAME}-1-2 sh -c "${BIN} hex1.mesh3d 1 x 2 y | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1-2.dump")
add_test(${PROJECT_NAME}-4-1 sh -c "${BIN} hex4.mesh3d 1 xyz | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/4-1.dump")
