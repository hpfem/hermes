include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/hermes_common/)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hermes1d/)

install(FILES __init__.py  DESTINATION ${PY_PATH})
install(FILES h1d_wrapper/__init__.py  DESTINATION ${PY_PATH}/h1d_wrapper)
install(FILES fekete/__init__.py  DESTINATION ${PY_PATH}/fekete)

file(GLOB H1D_WRAPPER_INC  "${CMAKE_CURRENT_SOURCE_DIR}/h1d_wrapper/*.pxd")
file(GLOB FEKETE_INC  "${CMAKE_CURRENT_SOURCE_DIR}/fekete/*.pxd")
install(FILES ${H1D_WRAPPER_INC}  DESTINATION ${PY_PATH}/h1d_wrapper)
install(FILES ${FEKETE_INC}  DESTINATION ${PY_PATH}/fekete)
