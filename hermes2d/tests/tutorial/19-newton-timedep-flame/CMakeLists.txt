if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(tutorial-19-newton-timedep-flame)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-19-newton-timedep-flame ${BIN})
set_tests_properties(tutorial-19-newton-timedep-flame PROPERTIES LABELS slow)
