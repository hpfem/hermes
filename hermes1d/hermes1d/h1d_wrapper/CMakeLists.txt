SET(CYTHON_INCLUDE_DIRECTORIES . -I../../../hermes_common/)
cython_add_module_pyx(h1d_wrapper)
add_python_library(h1d_wrapper h1d_wrapper.cpp)
target_link_libraries(h1d_wrapper ${HERMES_BIN})
