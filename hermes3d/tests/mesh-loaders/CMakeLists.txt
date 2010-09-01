project(mesh-loader)

if(H3D_REAL)

include(CMake.vars OPTIONAL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

add_executable(${PROJECT_NAME}
	main.cpp
	${HERMES_COMMON_DIR}/trace.cpp
)

target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES})

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

#
# Mesh3D
#

# cubes
if(WITH_HEX)
	add_test(mesh-loader-m3d-hex-1 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex1.dump")
	add_test(mesh-loader-m3d-hex-2 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex2.dump")
	add_test(mesh-loader-m3d-hex-3 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex4.dump")
	add_test(mesh-loader-m3d-hex-4 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex8.dump")
	
	add_test(mesh-loader-m3d-hex-5 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/missing_boundary.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/missing_boundary.dump")
	add_test(mesh-loader-m3d-hex-6 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/invalid_vertex_indices.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/invalid_vertex_indices.dump")
	add_test(mesh-loader-m3d-hex-7 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/missing_outer_face_boundary.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/missing_outer_face_boundary.dump")
endif(WITH_HEX)

# tetras
if(WITH_TETRA)
	add_test(mesh-loader-m3d-tet-1 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/tetra1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra1.dump")
	add_test(mesh-loader-m3d-tet-2 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/tetra2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra2.dump")
	add_test(mesh-loader-m3d-tet-3 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/tetra4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra4.dump")
	add_test(mesh-loader-m3d-tet-4 sh -c "${BIN} m3d ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra8.dump")
endif(WITH_TETRA)

#
# HDF5
#

if(WITH_HDF5)
# cubes
if(WITH_HEX)
	add_test(mesh-loader-hdf5-hex-1 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/hex1.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex1.dump")
	add_test(mesh-loader-hdf5-hex-2 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/hex2.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex2.dump")
	add_test(mesh-loader-hdf5-hex-3 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/hex4.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex4.dump")
	add_test(mesh-loader-hdf5-hex-4 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/hex8.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex8.dump")
endif(WITH_HEX)

# tetras
if(WITH_TETRA)
	add_test(mesh-loader-hdf5-tet-1 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/tetra1.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra1.dump")
	add_test(mesh-loader-hdf5-tet-2 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/tetra2.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra2.dump")
	add_test(mesh-loader-hdf5-tet-3 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/tetra4.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra4.dump")
	add_test(mesh-loader-hdf5-tet-4 sh -c "${BIN} hdf5 ${MESHES_DIR}/hdf5/tetra8.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/tetra8.dump")
endif(WITH_TETRA)

endif(WITH_HDF5)

#
# exodusII
#

if(WITH_EXODUSII)
# cubes
if(WITH_HEX)
	add_test(mesh-loader-exo2-hex-1 sh -c "${BIN} exoii ${MESHES_DIR}/exodusII/hex1.e | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex1-exoii.dump")
	add_test(mesh-loader-exo2-hex-2 sh -c "${BIN} exoii ${MESHES_DIR}/exodusII/hex27.e | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex27-exoii.dump")
	add_test(mesh-loader-exo2-hex-3 sh -c "${BIN} exoii ${MESHES_DIR}/exodusII/hex27-w-mat.e | diff - ${CMAKE_CURRENT_SOURCE_DIR}/dump/hex27-w-mat-exoii.dump")
endif(WITH_HEX)

endif(WITH_EXODUSII)

endif(H3D_REAL)
