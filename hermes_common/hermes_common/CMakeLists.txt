# TODO: Complex version?
if (NOT HERMES_COMMON_REAL)
    return()
endif (NOT HERMES_COMMON_REAL)

include_directories(${PROJECT_BINARY_DIR})

cython_add_module(cython_utils)
cython_add_module(numpy_utils)
cython_add_module(matrix)

# TODO: We currently allow only one version of matrix.so - if only the release
# version of hermes_common is being built, matrix.so will be linked with that
# version, otherwise it will be always linked with the debug version.
if (HERMES_COMMON_DEBUG)
    target_link_libraries(matrix ${HERMES_COMMON_LIB_REAL_DEBUG})
else (HERMES_COMMON_DEBUG)
    target_link_libraries(matrix ${HERMES_COMMON_LIB_REAL_RELEASE})
endif (HERMES_COMMON_DEBUG)

# FIXME: Seems that 'cython_utils' and 'matrix' aren't used outside H1D - verify
# and possibly make a switch here.
# TODO: Change 'hermes_common' to 'hermes/hermes_common'
set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes_common/)

install(TARGETS numpy_utils cython_utils matrix
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )

# This is needed, so that cmake doesn't install the symbolic link instead
# FIXME: Is this still true?
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
    
file(GLOB PY_INC  "${CMAKE_CURRENT_SOURCE_DIR}/*.pxd")
file(GLOB MAT_INC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/*")
install(FILES ${PY_INC}   DESTINATION ${PY_PATH})
install(FILES ${MAT_INC}  DESTINATION ${PY_PATH})
