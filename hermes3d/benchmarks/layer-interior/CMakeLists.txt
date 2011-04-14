project(h3d-layer-interior)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

