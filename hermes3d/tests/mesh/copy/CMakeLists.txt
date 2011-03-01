project(test-mesh-copy)

include(${hermes3d_SOURCE_DIR}/CMake.common)

add_executable(${PROJECT_NAME}	main.cpp)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

# copy
add_test(${PROJECT_NAME}-1-1 "${BIN}" copy hex1.mesh3d 1 x)
add_test(${PROJECT_NAME}-1-2 "${BIN}" copy hex1.mesh3d 1 y)
add_test(${PROJECT_NAME}-1-3 "${BIN}" copy hex1.mesh3d 1 z)

# copy base
add_test(${PROJECT_NAME}-base-1-1 "${BIN}" copy_base hex1.mesh3d)
add_test(${PROJECT_NAME}-base-1-2 "${BIN}" copy_base hex1.mesh3d 1 x)
add_test(${PROJECT_NAME}-base-1-3 "${BIN}" copy_base hex1.mesh3d 1 y)
add_test(${PROJECT_NAME}-base-1-4 "${BIN}" copy_base hex1.mesh3d 1 z)

add_test(${PROJECT_NAME}-base-2-1 "${BIN}" copy_base hex2.mesh3d)
add_test(${PROJECT_NAME}-base-2-2 "${BIN}" copy_base hex2.mesh3d 1 x)
add_test(${PROJECT_NAME}-base-2-3 "${BIN}" copy_base hex2.mesh3d 1 y)
add_test(${PROJECT_NAME}-base-2-4 "${BIN}" copy_base hex2.mesh3d 1 z)

add_test(${PROJECT_NAME}-base-4-1 "${BIN}" copy_base hex4.mesh3d)
add_test(${PROJECT_NAME}-base-4-2 "${BIN}" copy_base hex4.mesh3d 1 x)
add_test(${PROJECT_NAME}-base-4-3 "${BIN}" copy_base hex4.mesh3d 1 y)
add_test(${PROJECT_NAME}-base-4-4 "${BIN}" copy_base hex4.mesh3d 1 z)

add_test(${PROJECT_NAME}-base-8-1 "${BIN}" copy_base hex8.mesh3d)
add_test(${PROJECT_NAME}-base-8-2 "${BIN}" copy_base hex8.mesh3d 1 x)
add_test(${PROJECT_NAME}-base-8-3 "${BIN}" copy_base hex8.mesh3d 1 y)
add_test(${PROJECT_NAME}-base-8-4 "${BIN}" copy_base hex8.mesh3d 1 z)

add_test(${PROJECT_NAME}-base-27-1 "${BIN}" copy_base hex27.mesh3d)
add_test(${PROJECT_NAME}-base-27-2 "${BIN}" copy_base hex27.mesh3d 1 x)
add_test(${PROJECT_NAME}-base-27-3 "${BIN}" copy_base hex27.mesh3d 1 y)
add_test(${PROJECT_NAME}-base-27-4 "${BIN}" copy_base hex27.mesh3d 1 z)
