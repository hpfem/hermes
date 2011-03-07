project(test-h3d-mesh-loader)

include(${hermes3d_SOURCE_DIR}/CMake.common)

if(H3D_REAL)

add_executable(${PROJECT_NAME}	main.cpp)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_REAL})

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

#
# Mesh3D
#

# cubes
if(WITH_HEX)
  add_test(${PROJECT_NAME}-m3d-hex-1 "${BIN}" m3d hex1.mesh3d)
  add_test(${PROJECT_NAME}-m3d-hex-2 "${BIN}" m3d hex2.mesh3d)
  add_test(${PROJECT_NAME}-m3d-hex-3 "${BIN}" m3d hex4.mesh3d)
  add_test(${PROJECT_NAME}-m3d-hex-4 "${BIN}" m3d hex8.mesh3d)
  add_test(${PROJECT_NAME}-m3d-hex-6 "${BIN}" m3d invalid_vertex_indices.mesh3d supp_to_fail)
  if(HERMES_COMMON_CHECK_BOUNDARY_CONDITIONS)
    add_test(${PROJECT_NAME}-m3d-hex-5 "${BIN}" m3d missing_boundary.mesh3d supp_to_fail)
    add_test(${PROJECT_NAME}-m3d-hex-7 "${BIN}" m3d missing_outer_face_boundary.mesh3d supp_to_fail)
  endif(HERMES_COMMON_CHECK_BOUNDARY_CONDITIONS)
endif(WITH_HEX)

# tetras
if(WITH_TETRA)
  add_test(${PROJECT_NAME}-m3d-tet-1 "${BIN}" m3d tetra1.mesh3d)
  add_test(${PROJECT_NAME}-m3d-tet-2 "${BIN}" m3d tetra2.mesh3d)
  add_test(${PROJECT_NAME}-m3d-tet-3 "${BIN}" m3d tetra4.mesh3d)
  add_test(${PROJECT_NAME}-m3d-tet-4 "${BIN}" m3d tetra8.mesh3d)
endif(WITH_TETRA)

#
# HDF5
#

if(WITH_HDF5)
# cubes
if(WITH_HEX)
  add_test(${PROJECT_NAME}-hdf5-hex-1 "${BIN}" hdf5 hex1.h5)
  add_test(${PROJECT_NAME}-hdf5-hex-2 "${BIN}" hdf5 hex2.h5)
  add_test(${PROJECT_NAME}-hdf5-hex-3 "${BIN}" hdf5 hex4.h5)
  add_test(${PROJECT_NAME}-hdf5-hex-4 "${BIN}" hdf5 hex8.h5)
endif(WITH_HEX)

# tetras
if(WITH_TETRA)
  add_test(${PROJECT_NAME}-hdf5-tet-1 "${BIN}" hdf5 tetra1.h5)
  add_test(${PROJECT_NAME}-hdf5-tet-2 "${BIN}" hdf5 tetra2.h5)
  add_test(${PROJECT_NAME}-hdf5-tet-3 "${BIN}" hdf5 tetra4.h5)
  add_test(${PROJECT_NAME}-hdf5-tet-4 "${BIN}" hdf5 tetra8.h5)
endif(WITH_TETRA)

endif(WITH_HDF5)

#
# exodusII
#

if(WITH_EXODUSII)
# cubes
if(WITH_HEX)
  add_test(${PROJECT_NAME}-exo2-hex-1 "${BIN}" exoii hex1.e)
  add_test(${PROJECT_NAME}-exo2-hex-2 "${BIN}" exoii hex27.e)
  add_test(${PROJECT_NAME}-exo2-hex-3 "${BIN}" exoii hex27-w-mat.e)
endif(WITH_HEX)

endif(WITH_EXODUSII)

endif(H3D_REAL)
