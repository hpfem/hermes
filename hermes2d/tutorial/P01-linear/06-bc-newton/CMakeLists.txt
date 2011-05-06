project(P01-06-bc-newton)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)