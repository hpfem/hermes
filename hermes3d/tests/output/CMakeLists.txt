project(test-output)

if(H3D_REAL)

set(BIN_GMSH ${PROJECT_NAME}-gmsh)
set(BIN_VTK  ${PROJECT_NAME}-vtk)

include (${hermes3d_SOURCE_DIR}/CMake.common)

# GMSH ####

add_executable(${BIN_GMSH}
	main.cpp
)

set_common_target_properties(${BIN_GMSH})

get_target_property(CPFL ${BIN_GMSH} COMPILE_FLAGS)
set_target_properties(${BIN_GMSH}
    PROPERTIES
	COMPILE_FLAGS "${CPFL} -DGMSH")


# VTK ####

add_executable(${BIN_VTK}
	main.cpp
)

set_common_target_properties(${BIN_VTK})

get_target_property(CPFL ${BIN_VTK} COMPILE_FLAGS)
set_target_properties(${BIN_VTK}
    PROPERTIES
	COMPILE_FLAGS "${CPFL} -DVTK")

# Tests

# GMSH ####
set(BIN ${PROJECT_BINARY_DIR}/${BIN_GMSH})

if(WITH_HEX)
# solution
add_test(${PROJECT_NAME}-gmsh-sln-1 sh -c "${BIN} sln hex1.mesh3d")
add_test(${PROJECT_NAME}-gmsh-sln-2 sh -c "${BIN} sln hex2.mesh3d")
add_test(${PROJECT_NAME}-gmsh-sln-3 sh -c "${BIN} sln hex4.mesh3d")
add_test(${PROJECT_NAME}-gmsh-sln-4 sh -c "${BIN} sln hex8.mesh3d")
add_test(${PROJECT_NAME}-gmsh-sln-5 sh -c "${BIN} sln hex27.mesh3d")
add_test(${PROJECT_NAME}-gmsh-sln-6 sh -c "${BIN} sln fichera-corner.mesh3d")

# order
add_test(${PROJECT_NAME}-gmsh-ord-1 sh -c "${BIN} ord hex1.mesh3d")
add_test(${PROJECT_NAME}-gmsh-ord-2 sh -c "${BIN} ord hex2.mesh3d")
add_test(${PROJECT_NAME}-gmsh-ord-3 sh -c "${BIN} ord hex4.mesh3d")
add_test(${PROJECT_NAME}-gmsh-ord-4 sh -c "${BIN} ord hex8.mesh3d")
add_test(${PROJECT_NAME}-gmsh-ord-5 sh -c "${BIN} ord hex27.mesh3d")
add_test(${PROJECT_NAME}-gmsh-ord-6 sh -c "${BIN} ord fichera-corner.mesh3d")

# BC
add_test(${PROJECT_NAME}-gmsh-bc-5 sh -c "${BIN} bc hex27.mesh3d")
add_test(${PROJECT_NAME}-gmsh-bc-6 sh -c "${BIN} bc fichera-corner.mesh3d")

endif(WITH_HEX)


# VTK ####
set(BIN ${PROJECT_BINARY_DIR}/${BIN_VTK})

if(WITH_HEX)
# solution
add_test(${PROJECT_NAME}-vtk-sln-1 sh -c "${BIN} sln hex1.mesh3d")
add_test(${PROJECT_NAME}-vtk-sln-2 sh -c "${BIN} sln hex2.mesh3d")
add_test(${PROJECT_NAME}-vtk-sln-3 sh -c "${BIN} sln hex4.mesh3d")
add_test(${PROJECT_NAME}-vtk-sln-4 sh -c "${BIN} sln hex8.mesh3d")
add_test(${PROJECT_NAME}-vtk-sln-5 sh -c "${BIN} sln hex27.mesh3d")
add_test(${PROJECT_NAME}-vtk-sln-6 sh -c "${BIN} sln fichera-corner.mesh3d")

add_test(${PROJECT_NAME}-vtk-vec-sln-5 sh -c "${BIN} vec-sln hex27.mesh3d")

add_test(${PROJECT_NAME}-vtk-3sln-5 sh -c "${BIN} 3sln hex27.mesh3d")

# order
add_test(${PROJECT_NAME}-vtk-ord-1 sh -c "${BIN} ord hex1.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-2 sh -c "${BIN} ord hex2.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-3 sh -c "${BIN} ord hex4.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-4 sh -c "${BIN} ord hex8.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-5 sh -c "${BIN} ord hex27.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-6 sh -c "${BIN} ord fichera-corner.mesh3d")

# BC
add_test(${PROJECT_NAME}-vtk-bc-5 sh -c "${BIN} bc hex27.mesh3d")
add_test(${PROJECT_NAME}-vtk-bc-6 sh -c "${BIN} bc fichera-corner.mesh3d")

# mat
if(WITH_UMFPACK)	
add_test(${PROJECT_NAME}-vtk-mat-1 sh -c "${BIN} mat fichera-corner.mesh3d")
endif(WITH_UMFPACK)	

# multi-mesh
add_test(${PROJECT_NAME}-vtk-mm-1 sh -c "${BIN} mm hex1.mesh3d")

endif(WITH_HEX)

if(WITH_TETRA)
add_test(${PROJECT_NAME}-vtk-sln-tet-8 sh -c "${BIN} sln tetra8.mesh3d")
add_test(${PROJECT_NAME}-vtk-ord-tet-8 sh -c "${BIN} ord tetra8.mesh3d")
add_test(${PROJECT_NAME}-vtk-bc-tet-8 sh -c "${BIN} bc tetra8.mesh3d")
add_test(${PROJECT_NAME}-vtk-vec-sln-tet-8 sh -c "${BIN} vec-sln tetra8.mesh3d")
add_test(${PROJECT_NAME}-vtk-3sln-tet-8 sh -c "${BIN} 3sln tetra8.mesh3d")
endif(WITH_TETRA)

endif(H3D_REAL)
