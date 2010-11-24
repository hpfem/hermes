add_subdirectory(electrostatics)

set(PY_PATH lib/python/hermes2d/modules)
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
