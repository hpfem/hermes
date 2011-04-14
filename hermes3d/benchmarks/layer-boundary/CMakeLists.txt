project(h3d-layer-boundary)
add_executable(${PROJECT_NAME} main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# FIXME: This test was disabled on April 15, 2011
# since it took too long to finish. We need to 
# make it shorter.
#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

