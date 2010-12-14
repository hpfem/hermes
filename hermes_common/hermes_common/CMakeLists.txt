cython_add_module(cython_utils)
cython_add_module(numpy_utils)

cython_add_module_pyx(matrix)
add_python_library(matrix matrix.cpp)
target_link_libraries(matrix ${HERMES_BIN})

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/hermes_common/)
install(TARGETS numpy_utils
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
