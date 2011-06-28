if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)

project(P04-03-complex)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)