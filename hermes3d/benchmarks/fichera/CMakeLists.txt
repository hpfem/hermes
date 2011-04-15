project(fichera)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# FIXME - THIS TEST TOOK TOO LONG SO IT WAS 
# DISABLED ON APRIL 14, 2011.
#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

