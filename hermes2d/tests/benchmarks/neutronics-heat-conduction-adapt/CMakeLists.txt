if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-neutronics-heat-conduction-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-neutronics-heat-conduction-adapt ${BIN})
set_tests_properties(benchmark-neutronics-heat-conduction-adapt PROPERTIES LABELS slow)

