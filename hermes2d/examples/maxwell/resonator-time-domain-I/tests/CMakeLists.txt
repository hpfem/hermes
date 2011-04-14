
project(test-example-maxwell-resonator-time-domain-I)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_REAL})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-maxwell-resonator-time-domain-I ${BIN})

