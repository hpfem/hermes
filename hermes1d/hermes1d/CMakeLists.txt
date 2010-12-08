SET(CYTHON_INCLUDE_DIRECTORIES . -I${CMAKE_CURRENT_SOURCE_DIR}/hermes_common)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/../hermes_common/)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)
add_subdirectory(hermes_common)

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )
