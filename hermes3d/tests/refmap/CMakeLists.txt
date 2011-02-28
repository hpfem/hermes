project(test-ref-map)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME}-1 sh -c "${BIN} hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex1.dump")
add_test(${PROJECT_NAME}-2 sh -c "${BIN} hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex2.dump")
add_test(${PROJECT_NAME}-4 sh -c "${BIN} hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex4.dump")
add_test(${PROJECT_NAME}-8 sh -c "${BIN} hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex8.dump")
add_test(${PROJECT_NAME}-27 sh -c "${BIN} hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex27.dump")

