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
	ogprojection.cpp
	loader/exodusii.cpp
	loader/mesh3d.cpp
	loader/hdf5.cpp
	matrix.cpp
	norm.cpp
	output/gmsh.cpp
	output/vtk.cpp
	output/graph.cpp
	qsort.cpp
	quadcheb.cpp
	quadstd.cpp
	refdomain.cpp
	refmap.cpp
	shapefn.cpp
	shapeset/shapeset.cpp
	shapeset/lobatto.cpp
	shapeset/h1lobattotetra.cpp
	shapeset/h1lobattotetradx.cpp
	shapeset/h1lobattotetrady.cpp
	shapeset/h1lobattotetradz.cpp
	shapeset/h1lobattohex.cpp
	shapeset/hcurllobattohex.cpp
	shapeset/refmapss.cpp
	solution.cpp
	space/space.cpp
	space/h1.cpp
	space/hcurl.cpp
	transform.cpp
	traverse.cpp
	weakform.cpp
	solver/pardiso.cpp
	solver/petsc.cpp
	solver/umfpack_solver.cpp
	solver/epetra.cpp
	solver/amesos.cpp
	solver/aztecoo.cpp
	solver/nox.cpp
	solver/mumps.cpp
	solver/precond_ifpack.cpp
	solver/precond_ml.cpp
	adapt/adapt.cpp
	adapt/proj.cpp
	adapt/projipol.cpp
	adapt/h1proj.cpp
	adapt/h1projipol.cpp
	common_time_period.cpp
	${HERMES_COMMON_DIR}/callstack.cpp
	${HERMES_COMMON_DIR}/trace.cpp
	${HERMES_COMMON_DIR}/error.cpp
	${HERMES_COMMON_DIR}/utils.cpp
	${HERMES_COMMON_DIR}/timer.cpp
)

macro(LIBRARY_SETTINGS BIN)
	target_link_libraries(${BIN} ${JUDY_LIBRARY})
	target_link_libraries(${BIN} ${BLAS_LIBRARIES})
	if(NOT MSVC)
		IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(${BIN} "rt")
		ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	endif(NOT MSVC)

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
file(GLOB INC_OUTPUT    "${CMAKE_CURRENT_SOURCE_DIR}/output/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")

install(FILES ${INC}           DESTINATION include/)
install(FILES ${INC_COMMON}    DESTINATION include/common)
install(FILES ${INC_SPACE}     DESTINATION include/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/loader)
install(FILES ${INC_SOLVER}    DESTINATION include/solver)
install(FILES ${INC_OUTPUT}    DESTINATION include/output)
install(FILES ${INC_INTEGRALS} DESTINATION include/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/adapt)

# MSVC (Win) helper macros

# Makes Win32 path from Unix-style patch which is used by CMAKE. Used when a path is provided to an OS utility.
macro(MAKE_PATH PATH_OUT PATH_IN)
	if(WIN32)
		string(REPLACE "/" "\\" ${PATH_OUT} ${PATH_IN})
	else(WIN32)
		set(${PATH_OUT} ${PATH_IN})
	endif(WIN32)
endmacro(MAKE_PATH)

macro(BUILD_LIB_ALL LIB FLAGS)
	set_target_properties(${LIB} PROPERTIES DEBUG_OUTPUT_NAME "${LIB}")
	install(TARGETS ${LIB}
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION lib)
	MAKE_PATH(TARGET_DIR "${TARGET_ROOT}/bin")
	get_target_property(SOURCE_DEBUG_FILE ${LIB} LOCATION_Debug)
	MAKE_PATH(SOURCE_DEBUG_FILE ${SOURCE_DEBUG_FILE})
	get_target_property(SOURCE_RELEASE_FILE ${LIB} LOCATION_Release)
	MAKE_PATH(SOURCE_RELEASE_FILE ${SOURCE_RELEASE_FILE})
	add_custom_command(TARGET ${LIB}
		POST_BUILD
		COMMAND if exist ${SOURCE_DEBUG_FILE} copy /Y ${SOURCE_DEBUG_FILE} ${TARGET_DIR}
		COMMAND if exist ${SOURCE_RELEASE_FILE} copy /Y ${SOURCE_RELEASE_FILE} ${TARGET_DIR})
	unset(TARGET_DIR)
endmacro(BUILD_LIB_ALL)

# install library
if(H3D_REAL)
	if(MSVC)
		set(MSVC_DLL_DEFINES "-D_HERMESDLL")
		BUILD_LIB_ALL(${HERMES_REAL_BIN} "${MSVC_DLL_DEFINES}")	
	else(MSVC)
		install(TARGETS ${HERMES_REAL_BIN} LIBRARY DESTINATION lib)
	endif(MSVC)
endif(H3D_REAL)

if(H3D_COMPLEX)
	if(MSVC)
		set(MSVC_DLL_DEFINES "-D_HERMESDLL")
		BUILD_LIB_ALL(${HERMES_CPLX_BIN} "${MSVC_DLL_DEFINES} -DH2D_COMPLEX")	
	else(MSVC)
		install(TARGETS ${HERMES_CPLX_BIN} LIBRARY DESTINATION lib)
	endif(MSVC)
endif(H3D_COMPLEX)
