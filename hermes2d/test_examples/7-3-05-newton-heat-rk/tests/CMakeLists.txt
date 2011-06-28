if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(test-tutorial-P03-02-runge-kutta)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P03-02-cathedral-rk ${BIN})

