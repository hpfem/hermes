include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/../hermes_common/)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes1d/)

install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
file(GLOB PY_INC  "${CMAKE_CURRENT_SOURCE_DIR}/*.pxd")
install(FILES ${PY_INC}  DESTINATION ${PY_PATH})    
