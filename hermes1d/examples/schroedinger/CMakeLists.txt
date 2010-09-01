project(schroedinger)

find_package(PythonLibs REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hermes1d_SOURCE_DIR}/hermes1d/h1d_wrapper/)

add_executable(${PROJECT_NAME} main.cpp)
include(../CMake.common)
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
add_dependencies(${PROJECT_NAME} _hermes1d)
