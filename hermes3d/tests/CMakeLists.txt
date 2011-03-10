#
# CMake for tests
#
add_subdirectory(adapt)
add_subdirectory(calc)
add_subdirectory(hang-nodes)
add_subdirectory(mesh)
add_subdirectory(mesh-loaders)
add_subdirectory(orders)
add_subdirectory(output)
add_subdirectory(quadrature)
# This test needs some functions that are not implemented.
#add_subdirectory(refmap)
add_subdirectory(shapeset)
