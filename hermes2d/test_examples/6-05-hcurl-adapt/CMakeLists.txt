project(P04-04-hcurl)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})

if(H2D_WITH_TESTS)
  add_subdirectory(tests)
endif(H2D_WITH_TESTS)
