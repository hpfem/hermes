project(system_sin_eigen)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

if(WITH_TESTS)
  if(WITH_PYTHON)
#    FIXME: not adapted to the latest API
#    add_subdirectory(tests)
  endif(WITH_PYTHON)
endif(WITH_TESTS)

