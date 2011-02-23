IF(H2D_COMPLEX)

project(P03-08-gross-pitaevski)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)

ENDIF(H2D_COMPLEX)
