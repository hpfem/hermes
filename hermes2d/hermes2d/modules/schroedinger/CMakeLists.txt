set(CYTHON_INCLUDE_DIRECTORIES . -I ${PROJECT_BINARY_DIR}/hermes2d/hermes_common)
CYTHON_ADD_MODULE_PYX(schroedinger)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/modules/schroedinger)
add_python_library(schroedinger schroedinger.cpp)
target_link_libraries(schroedinger hermes_schroedinger)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/modules/schroedinger)
install(TARGETS schroedinger
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
