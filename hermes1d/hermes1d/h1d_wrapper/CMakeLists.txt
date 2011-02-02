cython_add_module_pyx(h1d_wrapper)
add_python_library(h1d_wrapper h1d_wrapper.cpp)
target_link_libraries(h1d_wrapper ${HERMES1D})