if(NOT H3D_REAL)
    return()
endif(NOT H3D_REAL)

project(test-example-elasticity-cubit-hex)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME} ${BIN})
