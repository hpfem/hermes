if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(tutorial-23-heat-and-moisture-adapt-long)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-23-heat-and-moisture-adapt-long ${BIN})
set_tests_properties(tutorial-23-heat-and-moisture-adapt-long PROPERTIES LABELS slow)
