#
# CMake for examples
#

# HEX
if(WITH_HEX)
	if(H3D_REAL)
		add_subdirectory(bridge)
		if(WITH_EXODUSII)
			add_subdirectory(elasticity-cubit-hex)
		endif(WITH_EXODUSII)
		add_subdirectory(heat-conduction)
		add_subdirectory(heat-cubit)
		add_subdirectory(schroedinger-eigen)
		add_subdirectory(hartree-fock-neweigen)
		add_subdirectory(schroedinger-eigen-refine-towards-origin)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		add_subdirectory(rough)
		add_subdirectory(maxwell-ring)
	endif(H3D_COMPLEX)
endif(WITH_HEX)

if(WITH_TETRA)
		add_subdirectory(elasticity-cubit-tet)
endif(WITH_TETRA)
