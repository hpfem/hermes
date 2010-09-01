project(hermes3d)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})		# for config.h

# source files for the hermes3d library
set(SRC
	hermes3d.cpp
	filter.cpp
	forms.cpp
	function.cpp
	mesh.cpp
	discrete_problem.cpp
	linear_problem.cpp
	loader/exodusii.cpp
	loader/mesh3d.cpp
	loader/hdf5.cpp
	matrix.cpp
	norm.cpp
	output/gmsh.cpp
	output/vtk.cpp
	output/graph.cpp
	precond/ifpack.cpp
	precond/ml.cpp
	qsort.cpp
	quadcheb.cpp
	quadstd.cpp
	refdomain.cpp
	refmap.cpp
	shapefn.cpp
	shapeset.cpp
	shapeset/lobatto.cpp
	shapeset/h1lobattotetra.cpp
	shapeset/h1lobattotetradx.cpp
	shapeset/h1lobattotetrady.cpp
	shapeset/h1lobattotetradz.cpp
	shapeset/h1lobattohex.cpp
	shapeset/hcurllobattohex.cpp
	shapeset/refmapss.cpp
	solution.cpp
	space.cpp
	space/h1.cpp
	space/hcurl.cpp
	transform.cpp
	traverse.cpp
	weakform.cpp
	solver/pardiso.cpp
	solver/petsc.cpp
	solver/umfpack.cpp
	solver/epetra.cpp
	solver/amesos.cpp
	solver/aztecoo.cpp
	solver/nox.cpp
	solver/mumps.cpp
	adapt/h1.cpp
	adapt/proj.cpp
	adapt/projipol.cpp
	adapt/h1proj.cpp
	adapt/h1projipol.cpp
	${HERMES_COMMON_DIR}/callstack.cpp
	${HERMES_COMMON_DIR}/trace.cpp
	${HERMES_COMMON_DIR}/error.cpp
	${HERMES_COMMON_DIR}/utils.cpp
	${HERMES_COMMON_DIR}/timer.cpp
)

macro(LIBRARY_SETTINGS BIN)
	target_link_libraries(${BIN} ${JUDY_LIBRARY})
	target_link_libraries(${BIN} ${BLAS_LIBRARIES})
	if(WITH_PETSC)
		include_directories(${PETSC_INCLUDE_DIR})
		target_link_libraries(${BIN} ${PETSC_LIBRARIES})
		target_link_libraries(${BIN} ${PETSC_LIBS})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_PETSC)

	if(WITH_UMFPACK)
		include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
		target_link_libraries(${BIN} ${UMFPACK_LIBRARY} ${AMD_LIBRARY})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_UMFPACK)

	if(WITH_PARDISO)
		target_link_libraries(${BIN} ${PARDISO_LIBRARY})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_PARDISO)

	if(WITH_MUMPS)
		include_directories(${MUMPS_INCLUDE_DIR})
		target_link_libraries(${BIN} ${MUMPS_LIBRARIES})
		target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_MUMPS)

	if(WITH_TRILINOS)
		include_directories(${TRILINOS_INCLUDE_DIR})
		target_link_libraries(${BIN} ${TRILINOS_LIBRARIES})
	endif(WITH_TRILINOS)

	if(WITH_HDF5)
	endif(WITH_HDF5)

	if(WITH_EXODUSII)
		include_directories(${EXODUSII_INCLUDE_DIR})
		target_link_libraries(${BIN} ${EXODUSII_LIBRARIES})
	endif(WITH_EXODUSII)

	target_link_libraries(${BIN} ${ADDITIONAL_LIBS})
endmacro(LIBRARY_SETTINGS)

# real version
if(H3D_REAL)
	add_library(${HERMES_REAL_BIN} SHARED ${SRC})
	LIBRARY_SETTINGS(${HERMES_REAL_BIN})
endif(H3D_REAL)

# complex version
if(H3D_COMPLEX)
	add_library(${HERMES_CPLX_BIN} SHARED ${SRC})
	set_target_properties(${HERMES_CPLX_BIN} PROPERTIES COMPILE_FLAGS
-DH3D_COMPLEX)
	LIBRARY_SETTINGS(${HERMES_CPLX_BIN})
endif(H3D_COMPLEX)

#
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/h3dconfig.h
)


#
# INSTALL
#

# install header files
file(GLOB INC           "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB INC_COMMON    "${HERMES_COMMON_DIR}/*.h")
file(GLOB INC_SPACE     "${CMAKE_CURRENT_SOURCE_DIR}/space/*.h")
file(GLOB INC_SHAPESET  "${CMAKE_CURRENT_SOURCE_DIR}/shapeset/*.h")
file(GLOB INC_LOADER    "${CMAKE_CURRENT_SOURCE_DIR}/loader/*.h")
file(GLOB INC_SOLVER    "${CMAKE_CURRENT_SOURCE_DIR}/solver/*.h")
file(GLOB INC_PRECOND   "${CMAKE_CURRENT_SOURCE_DIR}/precond/*.h")
file(GLOB INC_OUTPUT    "${CMAKE_CURRENT_SOURCE_DIR}/output/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")

install(FILES ${INC}           DESTINATION include/)
install(FILES ${INC_COMMON}    DESTINATION include/common)
install(FILES ${INC_SPACE}     DESTINATION include/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/loader)
install(FILES ${INC_SOLVER}    DESTINATION include/solver)
install(FILES ${INC_PRECOND}   DESTINATION include/precond)
install(FILES ${INC_OUTPUT}    DESTINATION include/output)
install(FILES ${INC_INTEGRALS} DESTINATION include/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/adapt)

# install library
if(H3D_REAL)
	install(TARGETS ${HERMES_REAL_BIN} LIBRARY DESTINATION lib)
endif(H3D_REAL)

if(H3D_COMPLEX)
	install(TARGETS ${HERMES_CPLX_BIN} LIBRARY DESTINATION lib)
endif(H3D_COMPLEX)
