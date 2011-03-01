project(test-projection)
add_executable(${PROJECT_NAME}	main.cpp)

include(${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

#
# Tests
#

enable_testing()

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-h1-1 ${BIN} h1 hex1.mesh3d)
add_test(${PROJECT_NAME}-h1-2 ${BIN} h1 hex2.mesh3d)
add_test(${PROJECT_NAME}-h1-3 ${BIN} h1 hex4.mesh3d)
add_test(${PROJECT_NAME}-h1-4 ${BIN} h1 hex8.mesh3d)
add_test(${PROJECT_NAME}-h1-5 ${BIN} h1 hex27.mesh3d)

add_test(${PROJECT_NAME}-h1i-1 ${BIN} h1-ipol hex1.mesh3d)
add_test(${PROJECT_NAME}-h1i-2 ${BIN} h1-ipol hex2.mesh3d)
add_test(${PROJECT_NAME}-h1i-3 ${BIN} h1-ipol hex4.mesh3d)
add_test(${PROJECT_NAME}-h1i-4 ${BIN} h1-ipol hex8.mesh3d)
add_test(${PROJECT_NAME}-h1i-5 ${BIN} h1-ipol hex27.mesh3d)

set_tests_properties(${PROJECT_NAME}-h1-5 PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-h1i-4 PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-h1i-5 PROPERTIES LABELS slow)
