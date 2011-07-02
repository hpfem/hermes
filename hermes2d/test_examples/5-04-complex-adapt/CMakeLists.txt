project(P04-03-complex)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})

if(H2D_WITH_TESTS)
  add_subdirectory(tests)
endif(H2D_WITH_TESTS)
