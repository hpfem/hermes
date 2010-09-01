SET(CYTHON_INCLUDE_DIRECTORIES . -I../../../hermes_common/)
cython_add_module(h1d_wrapper)
target_link_libraries(h1d_wrapper hermes_common ${HERMES_BIN})
set_property(TARGET h1d_wrapper PROPERTY COMPILE_FLAGS "-g")

install(TARGETS h1d_wrapper
            RUNTIME DESTINATION lib/python/hermes1d/h1d_wrapper/
            ARCHIVE DESTINATION lib/python/hermes1d/h1d_wrapper/
            LIBRARY DESTINATION lib/python/hermes1d/h1d_wrapper/
        )
set_target_properties(h1d_wrapper PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes1d/h1d_wrapper/"
    )

install(FILES __init__.py
    DESTINATION lib/python/hermes1d/h1d_wrapper/
    )

install(FILES _hermes1d.pxd __init__.py
    DESTINATION lib/python/hermes1d/include/hermes1d/h1d_wrapper/
    )
