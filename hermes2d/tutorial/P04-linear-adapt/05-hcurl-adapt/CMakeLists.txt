IF(H2D_COMPLEX)

project(P04-05-hcurl-adapt)

add_executable(${PROJECT_NAME} main.cpp bessel.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)

ENDIF(H2D_COMPLEX)
