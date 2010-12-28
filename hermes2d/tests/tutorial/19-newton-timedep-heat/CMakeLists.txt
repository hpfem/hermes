project(tutorial-18-newton-timedep-heat)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-18-newton-timedep-heat ${BIN})
set_tests_properties(tutorial-18-newton-timedep-heat PROPERTIES LABELS slow)
