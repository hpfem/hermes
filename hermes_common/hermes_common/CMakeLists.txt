cython_add_module_pyx(matrix)
add_python_library(matrix matrix.cpp)
target_link_libraries(matrix ${HERMES_BIN})
