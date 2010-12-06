CYTHON_ADD_MODULE_PYX(basic)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/modules/basic)
add_python_library(basic basic.cpp)
target_link_libraries(basic hermes_basic)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d/modules/basic)
install(TARGETS basic
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
