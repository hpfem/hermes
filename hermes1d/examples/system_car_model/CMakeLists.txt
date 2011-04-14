project(system_car_model)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

# FIXME: test disabled because it was too slow. Revision is needed.
#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

