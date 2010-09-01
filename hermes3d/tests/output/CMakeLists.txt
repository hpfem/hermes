project(output)

if(H3D_REAL)

include(CMake.vars OPTIONAL)

set(BIN_GMSH ${PROJECT_NAME}-gmsh)
set(BIN_VTK  ${PROJECT_NAME}-vtk)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

# GMSH ####

add_executable(${BIN_GMSH}
	main.cpp
)

set_target_properties(${BIN_GMSH}
    PROPERTIES
	COMPILE_FLAGS -DGMSH)


target_link_libraries(${BIN_GMSH} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_GMSH} ${HERMES_REAL_BIN})

# VTK ####

add_executable(${BIN_VTK}
	main.cpp
)

set_target_properties(${BIN_VTK}
    PROPERTIES
	COMPILE_FLAGS -DVTK)

target_link_libraries(${BIN_VTK} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_VTK} ${HERMES_REAL_BIN})


configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

# GMSH ####
set(BIN ${PROJECT_BINARY_DIR}/${BIN_GMSH})

if(WITH_HEX)
# solution
add_test(${PROJECT_NAME}-gmsh-sln-1 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-1")
add_test(${PROJECT_NAME}-gmsh-sln-2 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-2")
add_test(${PROJECT_NAME}-gmsh-sln-3 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-4")
add_test(${PROJECT_NAME}-gmsh-sln-4 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-8")
add_test(${PROJECT_NAME}-gmsh-sln-5 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-27")
add_test(${PROJECT_NAME}-gmsh-sln-6 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-fichera-corner")

# order
add_test(${PROJECT_NAME}-gmsh-ord-1 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-1")
add_test(${PROJECT_NAME}-gmsh-ord-2 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-2")
add_test(${PROJECT_NAME}-gmsh-ord-3 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-4")
add_test(${PROJECT_NAME}-gmsh-ord-4 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-8")
add_test(${PROJECT_NAME}-gmsh-ord-5 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-27")
add_test(${PROJECT_NAME}-gmsh-ord-6 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-fichera-corner")

# BC
add_test(${PROJECT_NAME}-gmsh-bc-5 sh -c "${BIN} bc ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-bc-27")
add_test(${PROJECT_NAME}-gmsh-bc-6 sh -c "${BIN} bc ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-bc-fichera-corner")

endif(WITH_HEX)


# VTK ####
set(BIN ${PROJECT_BINARY_DIR}/${BIN_VTK})

if(WITH_HEX)
# solution
add_test(${PROJECT_NAME}-vtk-sln-1 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-1")
add_test(${PROJECT_NAME}-vtk-sln-2 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-2")
add_test(${PROJECT_NAME}-vtk-sln-3 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-4")
add_test(${PROJECT_NAME}-vtk-sln-4 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-8")
add_test(${PROJECT_NAME}-vtk-sln-5 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-27")
add_test(${PROJECT_NAME}-vtk-sln-6 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-fichera-corner")

add_test(${PROJECT_NAME}-vtk-vec-sln-5 sh -c "${BIN} vec-sln ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-vec-sln-27")

add_test(${PROJECT_NAME}-vtk-3sln-5 sh -c "${BIN} 3sln ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-3sln-27")

# order
add_test(${PROJECT_NAME}-vtk-ord-1 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-1")
add_test(${PROJECT_NAME}-vtk-ord-2 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-2")
add_test(${PROJECT_NAME}-vtk-ord-3 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-4")
add_test(${PROJECT_NAME}-vtk-ord-4 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-8")
add_test(${PROJECT_NAME}-vtk-ord-5 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-27")
add_test(${PROJECT_NAME}-vtk-ord-6 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-fichera-corner")

# BC
add_test(${PROJECT_NAME}-vtk-bc-5 sh -c "${BIN} bc ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-bc-27")
add_test(${PROJECT_NAME}-vtk-bc-6 sh -c "${BIN} bc ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-bc-fichera-corner")

# mat
if(WITH_UMFPACK)	
add_test(${PROJECT_NAME}-vtk-mat-1 sh -c "${BIN} mat ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-mat-fichera-corner")
endif(WITH_UMFPACK)	

# multi-mesh
add_test(${PROJECT_NAME}-vtk-mm-1 sh -c "${BIN} mm ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-mm-1")

endif(WITH_HEX)

if(WITH_TETRA)
add_test(${PROJECT_NAME}-vtk-sln-tet-8 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-sln-tet-8")
add_test(${PROJECT_NAME}-vtk-ord-tet-8 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-ord-tet-8")
add_test(${PROJECT_NAME}-vtk-bc-tet-8 sh -c "${BIN} bc ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-bc-tet-8")
add_test(${PROJECT_NAME}-vtk-vec-sln-tet-8 sh -c "${BIN} vec-sln ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-vec-sln-tet-8")
add_test(${PROJECT_NAME}-vtk-3sln-tet-8 sh -c "${BIN} 3sln ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/vtk-3sln-tet-8")
endif(WITH_TETRA)

endif(H3D_REAL)
