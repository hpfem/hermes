#
# CMake for examples
#

include (../tests/common.cmake)

# HEX
if(WITH_HEX)
    # H1
	if(H3D_REAL)
		add_subdirectory(fichera)
		add_subdirectory(layer-boundary)
		add_subdirectory(layer-interior)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		add_subdirectory(bessel)
	endif(H3D_COMPLEX)
endif(WITH_HEX)
