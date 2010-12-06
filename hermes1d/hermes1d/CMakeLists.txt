set(HERMES_COMMON_PATH ${PROJECT_SOURCE_DIR}/../hermes_common)
execute_process(
    COMMAND cp -r ${HERMES_COMMON_PATH}/hermes_common ${CMAKE_CURRENT_SOURCE_DIR}/
    )
SET(CYTHON_INCLUDE_DIRECTORIES . -I${CMAKE_CURRENT_SOURCE_DIR}/hermes_common)

include_directories(${PROJECT_BINARY_DIR}/src)

cython_add_module(cython_utils)
cython_add_module(numpy_utils)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)
add_subdirectory(hermes_common)

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )
