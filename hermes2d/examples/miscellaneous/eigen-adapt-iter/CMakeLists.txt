
project(eigen-adapt-iter)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

add_subdirectory(tests)