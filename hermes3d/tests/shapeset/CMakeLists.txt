#
# Testing shapesets
#
if(H3D_REAL)
	# test cases with H3D_REAL version od Hermes3D
	add_subdirectory(hex-h1)
	#add_subdirectory(hex-h1-lobatto-ced)
	add_subdirectory(hex-hcurl)
endif(H3D_REAL)
