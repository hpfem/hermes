#
# CMake for test the calculations
#

# TETRA
if(WITH_TETRA)
	if(H3D_REAL)
		add_subdirectory(tetra-h1)
		add_subdirectory(tetra-h1-dirichlet)
		add_subdirectory(tetra-h1-neumann)
		add_subdirectory(tetra-h1-newton)
	endif(H3D_REAL)
endif(WITH_TETRA)

# HEX
if(WITH_HEX)
    # H1
	if(H3D_REAL)
		add_subdirectory(hex-h1)
		add_subdirectory(hex-h1-1face)
		add_subdirectory(hex-h1-dirichlet)
		add_subdirectory(hex-h1-neumann)
		add_subdirectory(hex-h1-newton)
		add_subdirectory(hex-h1-unsym)
		# systems of equations
		add_subdirectory(hex-h1-sys)
		add_subdirectory(hex-h1-sys-dirichlet)
		add_subdirectory(hex-h1-sys-coupled)
		# non linear
		if(HAVE_NOX)
			add_subdirectory(hex-h1-nonlin)
			add_subdirectory(hex-h1-matfree)
		endif(HAVE_NOX)
	endif(H3D_REAL)
	# Hcurl
	if(H3D_REAL)
		add_subdirectory(hex-hcurl)		
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		add_subdirectory(hex-hcurl-cplx)
		add_subdirectory(hex-hcurl-imped)
	endif(H3D_COMPLEX)
endif(WITH_HEX)
