project(test-tutorial-P03-03-newton-heat-ie)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P03-03-newton-heat-ie ${BIN})
set_tests_properties(test-tutorial-P03-03-newton-heat-ie PROPERTIES LABELS slow)
