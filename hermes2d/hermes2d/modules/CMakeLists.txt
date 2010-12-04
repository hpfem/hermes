add_subdirectory(basic)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/modules)
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
