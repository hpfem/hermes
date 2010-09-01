set(CYTHON_INCLUDE_DIRECTORIES ../..)

include_directories(..)

macro(ADD_EXAMPLE name)
    cython_add_module(${name})
    target_link_libraries(${name} ${HERMES_REAL_BIN})
    if(WITH_PYTHON_UMFPACK)
        find_package(UMFPACK REQUIRED)
        include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
        target_link_libraries(${name} umfpack amd blas)
    endif(WITH_PYTHON_UMFPACK)
    set_target_properties(${name} PROPERTIES
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes2d/examples"
        )
endmacro(ADD_EXAMPLE)

add_example(c03)
add_example(c04)
add_example(c05)
add_example(c06)
add_example(c07)
add_example(c08)
add_example(c09)
add_example(c10)
add_example(c11)
add_example(c12)
add_example(c22)

install(TARGETS c03 c04 c05 c06 c07 c08 c09 c10 c11 c12 c22
            RUNTIME DESTINATION lib/python/hermes2d/examples
            ARCHIVE DESTINATION lib/python/hermes2d/examples
            LIBRARY DESTINATION lib/python/hermes2d/examples
        )

install(FILES __init__.py domain.mesh sample.mesh cylinder4.mesh 07.mesh
    12.mesh bracket.mesh cathedral.mesh motor.mesh
    DESTINATION lib/python/hermes2d/examples
    )

