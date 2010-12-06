project(cppschroedinger)

find_package(PythonLibs REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hermes1d_SOURCE_DIR}/hermes1d/h1d_wrapper/)

set(cppschroedinger ${PROJECT_NAME})
add_executable(${cppschroedinger} main.cpp forms.cpp)
include(../CMake.common)
target_link_libraries(${cppschroedinger} ${PYTHON_LIBRARIES})
add_dependencies(${cppschroedinger} _hermes1d)

SET(CYTHON_INCLUDE_DIRECTORIES . -I../../hermes_common/ -I../../)
cython_add_module(_forms forms.cpp)
target_link_libraries(_forms ${HERMES_BIN})
