include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/../hermes_common/)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )
