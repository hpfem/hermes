cython_add_module(cython_utils)
cython_add_module(numpy_utils)

cython_add_module_pyx(matrix)
add_python_library(matrix matrix.cpp)
target_link_libraries(matrix ${HERMES})

# FIXME: this is dimension dependent:
set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/hermes_common/)
install(TARGETS numpy_utils
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
# This is needed, so that cmake doesn't install the symbolic link instead:
install(FILES ${PROJECT_SOURCE_DIR}/../hermes_common/hermes_common/__init__.py
    DESTINATION ${PY_PATH}
    )
