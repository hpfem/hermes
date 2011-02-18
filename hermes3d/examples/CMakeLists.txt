#
# CMake for examples
#

# HEX
if(WITH_HEX)
    # H1
	if(H3D_REAL)
		add_subdirectory(elasticity-cubit-hex)
		add_subdirectory(heat-conduction)
		add_subdirectory(heat-cubit)
		add_subdirectory(schroedinger-eigen)
		add_subdirectory(molecule3d-eigen)
		add_subdirectory(hartree-fock-eigen)
		add_subdirectory(hartree-fock-neweigen)
		#add_subdirectory(hartree-fock-eigen-new)
		#add_subdirectory(schroedinger-eigen-cusp)
		add_subdirectory(schroedinger-eigen-refine-towards-origin)
		#add_subdirectory(poisson-equation-box)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		add_subdirectory(rough)
		add_subdirectory(maxwell-ring)
	endif(H3D_COMPLEX)
endif(WITH_HEX)

if(WITH_TETRA)
		add_subdirectory(elasticity-cubit-tet)
endif(WITH_TETRA)
