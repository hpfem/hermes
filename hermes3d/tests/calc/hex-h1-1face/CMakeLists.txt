project(test-hex-h1-1face)

if(DEV_TESTS)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

enable_testing()
add_test(${PROJECT_NAME}-1 ${PROJECT_NAME})

endif(DEV_TESTS)
