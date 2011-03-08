add_subdirectory(quadrature)
add_subdirectory(bubbles)
add_subdirectory(mesh)
add_subdirectory(adaptivity)
add_subdirectory(shapeset)
add_subdirectory(integrals)
add_subdirectory(rcp)
add_subdirectory(python)

# Additional definitions for tests.
add_definitions(-DHERMES_REPORT_ALL -DH2D_TEST)
