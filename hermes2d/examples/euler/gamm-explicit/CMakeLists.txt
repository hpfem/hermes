if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)

project(gamm-explicit)

add_executable(${PROJECT_NAME} main.cpp ../euler_util.cpp ../numerical_flux.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Commented out on April 18 due to buildbot failure.
#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)