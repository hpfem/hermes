project(cppdirac)

include_directories(${hermes1d_SOURCE_DIR}/src)
add_executable(${PROJECT_NAME} main.cpp forms.cpp)

SET(CYTHON_INCLUDE_DIRECTORIES -I${CMAKE_SOURCE_DIR}/hermes_common)
SET(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${hermes1d_SOURCE_DIR})

cython_add_module(forms_wrapper forms.cpp)
target_link_libraries(forms_wrapper ${HERMES1D})
