add_subdirectory(h1d_wrapper)
add_subdirectory(fekete)

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )
