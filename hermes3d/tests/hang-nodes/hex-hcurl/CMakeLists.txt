project(hnnd-hex-hcurl)

add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
set(NAME ${PROJECT_NAME})
include(../tests-hex.cmake)

set_tests_properties(hnnd-hex-hcurl-1-01 PROPERTIES LABELS slow)
set_tests_properties(hnnd-hex-hcurl-1-02 PROPERTIES LABELS slow)
