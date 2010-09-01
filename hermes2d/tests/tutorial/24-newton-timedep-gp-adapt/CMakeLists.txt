if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)
project(24-newton-timedep-gp-adapt)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-24-newton-timedep-gp-adapt ${BIN})
set_tests_properties(tutorial-24-newton-timedep-gp-adapt PROPERTIES LABELS slow)

