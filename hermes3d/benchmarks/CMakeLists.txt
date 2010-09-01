#
# CMake for examples
#

include (../tests/common.cmake)

# HEX
if(WITH_HEX)
    # H1
	if(H3D_REAL)
#		add_subdirectory(elastostatics)
		add_subdirectory(fichera)
#		add_subdirectory(heat-conduction)
#		add_subdirectory(heat-cubit)
		add_subdirectory(layer)
#		add_subdirectory(singpert-aniso)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		add_subdirectory(bessel)
	endif(H3D_COMPLEX)
endif(WITH_HEX)
