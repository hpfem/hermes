if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(tutorial-22-newton-timedep-heat-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-22-newton-timedep-heat-adapt ${BIN})
set_tests_properties(tutorial-22-newton-timedep-heat-adapt PROPERTIES LABELS slow)
