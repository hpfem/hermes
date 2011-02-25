include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/../hermes_common/)

#add_subdirectory(hermes_common)

add_python_library(hermes2d ../src/hermes2d.cpp)
include(PickRealOrCplxLibs)
PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON_LIB_REAL} hermes2d)
target_link_libraries(hermes2d ${HERMES2D_REAL})
set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes2d)
install(TARGETS hermes2d
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py plot.py
    DESTINATION ${PY_PATH}
    )

add_subdirectory(hermes_common)
