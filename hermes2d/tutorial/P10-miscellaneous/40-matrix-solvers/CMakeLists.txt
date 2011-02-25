
if (H2D_REAL)
  project(P10-40-matrix-solvers-real)
  add_executable(${PROJECT_NAME} main.cpp)
  include (${hermes2d_SOURCE_DIR}/CMake.common)
  set_common_target_properties(${PROJECT_NAME} ${HERMES2D_REAL})
endif(H2D_REAL)

if (H2D_COMPLEX)
  project(P10-40-matrix-solvers-cplx)
  add_executable(${PROJECT_NAME} main.cpp)
  include (${hermes2d_SOURCE_DIR}/CMake.common)
  set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
endif(H2D_COMPLEX)

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)