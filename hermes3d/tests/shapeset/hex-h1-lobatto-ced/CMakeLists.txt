project(test-hex-h1-lobatto-shapeset-ced)

if(WITH_HEX)
#
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN})
set_tests_properties(${PROJECT_NAME}-1 PROPERTIES LABELS slow)

endif(WITH_HEX)

