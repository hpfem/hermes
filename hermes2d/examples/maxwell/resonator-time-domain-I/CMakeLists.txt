
project(resonator-time-domain-I)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_REAL})

#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

