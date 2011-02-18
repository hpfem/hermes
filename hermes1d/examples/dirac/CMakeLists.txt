project(cppdirac)

find_package(PythonLibs REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hermes1d_SOURCE_DIR}/hermes1d/h1d_wrapper/)

add_executable(${PROJECT_NAME} main.cpp forms.cpp)



SET(CYTHON_INCLUDE_DIRECTORIES . -I../../hermes_common/ -I../../)
cython_add_module(forms_wrapper forms.cpp)
target_link_libraries(forms_wrapper ${HERMES1D})
