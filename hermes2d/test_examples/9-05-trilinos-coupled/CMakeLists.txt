if(WITH_TRILINOS)
  project(P09-05-trilinos-coupled)

  add_executable(${PROJECT_NAME} main.cpp)
  include (${hermes2d_SOURCE_DIR}/CMake.common)
  set_common_target_properties(${PROJECT_NAME})

  if(H2D_WITH_TESTS)
    add_subdirectory(tests)
  endif(H2D_WITH_TESTS)
endif(WITH_TRILINOS)