project(P03-02-runge-kutta)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

if(H2D_WITH_TESTS)
  add_subdirectory(tests)
endif(H2D_WITH_TESTS)