add_subdirectory(examples)
SET(CYTHON_INCLUDE_DIRECTORIES . -I../../hermes_common/)
cython_add_module(_hermes2d dummy_solver.cpp stdcython.c)
add_custom_target(convert_api ./convert_api.py
    DEPENDS _hermes2d.cpp
    )
add_dependencies(_hermes2d convert_api)
target_link_libraries(_hermes2d hermes_common ${HERMES2D_REAL})
set_property(TARGET _hermes2d PROPERTY COMPILE_FLAGS "-g -DH2D_REPORT_RUNTIME_CONTROL")
if(WITH_PYTHON_UMFPACK)
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
    target_link_libraries(_hermes2d umfpack amd blas)
    set(USE_UMFPACK 1)
else(WITH_PYTHON_UMFPACK)
    set(USE_UMFPACK 0)
endif(WITH_PYTHON_UMFPACK)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

cython_add_module(forms)
target_link_libraries(forms ${HERMES2D_REAL})
install(TARGETS forms
            RUNTIME DESTINATION lib/python/hermes2d
            ARCHIVE DESTINATION lib/python/hermes2d
            LIBRARY DESTINATION lib/python/hermes2d
        )

set_target_properties(forms PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes2d"
    )

install(TARGETS _hermes2d
            RUNTIME DESTINATION lib/python/hermes2d
            ARCHIVE DESTINATION lib/python/hermes2d
            LIBRARY DESTINATION lib/python/hermes2d
        )
set_target_properties(_hermes2d PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes2d"
    )

install(FILES __init__.py plot.py runtests.py mesh.py demos.py
    DESTINATION lib/python/hermes2d
    )

# tests
install(FILES tests/test_basics.py tests/test_examples.py tests/test_plots.py
    tests/test_mesh.py
    DESTINATION lib/python/hermes2d/tests
    )

install(FILES dummy_solver.h config.h stdcython.h
    DESTINATION lib/python/hermes2d/include
    )

install(FILES _hermes2d.pxd __init__.py
    DESTINATION lib/python/hermes2d/include/hermes2d
    )
