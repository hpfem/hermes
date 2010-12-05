include_directories(${PROJECT_BINARY_DIR}/src)

cython_add_module(cython_utils)
cython_add_module(numpy_utils)

add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )
