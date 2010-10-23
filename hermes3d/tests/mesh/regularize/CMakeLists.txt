project(regularize)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1 sh -c "${BIN} hex4.mesh3d 2 xyz | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/1.dump")
add_test(${PROJECT_NAME}-2 sh -c "${BIN} hex4.mesh3d 2 xyz 5 xyz | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/2.dump")
