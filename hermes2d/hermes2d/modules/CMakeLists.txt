# Disable basic, as it doesn't compile:
#add_subdirectory(basic)
add_subdirectory(schroedinger)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/modules)
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
