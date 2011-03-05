include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/../hermes_common/)  # TODO: probably not needed

# TODO: Decide whether to link to DEBUG or RELEASE version of Hermes ?

cython_add_module(hermes2d)

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
file(GLOB PY_INC  "${CMAKE_CURRENT_SOURCE_DIR}/*.pxd")
install(FILES ${PY_INC}  DESTINATION ${PY_PATH})    
