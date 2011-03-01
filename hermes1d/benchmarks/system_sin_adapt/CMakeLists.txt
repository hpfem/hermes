project(system_sin_adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)

