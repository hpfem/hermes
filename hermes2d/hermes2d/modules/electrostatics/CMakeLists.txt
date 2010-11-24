CYTHON_ADD_MODULE_PYX(electrostatics)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/modules/electrostatics)
add_python_library(electrostatics electrostatics.cpp)
target_link_libraries(electrostatics hermes_electrostatics)

set(PY_PATH lib/python/hermes2d/modules/electrostatics)
install(TARGETS electrostatics
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
