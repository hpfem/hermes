# Enabling multiprocessor build on MSVC
if(MSVC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif(MSVC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# source files for the hermes3d library
set(SRC
	hermes3d.cpp
	filter.cpp
	weakform/forms.cpp
	function.cpp
	mesh.cpp
	discrete_problem.cpp
	ogprojection.cpp
	loader/ctuReader.cpp
	loader/exodusii.cpp
	loader/mesh3d.cpp
	loader/hdf5.cpp
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
	weakform/weakform.cpp
        adapt/adapt.cpp
	adapt/proj.cpp
	adapt/projipol.cpp
	adapt/h1proj.cpp
	adapt/hcurlproj.cpp
	adapt/h1projipol.cpp
  ${HERMES_COMMON_DIR}/common.cpp
  ${HERMES_COMMON_DIR}/hermes_logging.cpp
  ${HERMES_COMMON_DIR}/common_time_period.cpp
  ${HERMES_COMMON_DIR}/callstack.cpp
  ${HERMES_COMMON_DIR}/error.cpp
  ${HERMES_COMMON_DIR}/utils.cpp
  ${HERMES_COMMON_DIR}/matrix.cpp
  ${HERMES_COMMON_DIR}/tables.cpp
  ${HERMES_COMMON_DIR}/trace.cpp
  ${HERMES_COMMON_DIR}/third_party_codes/trilinos-teuchos/Teuchos_stacktrace.cpp 
  ${HERMES_COMMON_DIR}/solver/nox.cpp 
  ${HERMES_COMMON_DIR}/solver/epetra.cpp 
  ${HERMES_COMMON_DIR}/solver/aztecoo.cpp 
  ${HERMES_COMMON_DIR}/solver/amesos.cpp
  ${HERMES_COMMON_DIR}/solver/mumps.cpp 
  ${HERMES_COMMON_DIR}/solver/pardiso.cpp 
  ${HERMES_COMMON_DIR}/solver/petsc.cpp 
  ${HERMES_COMMON_DIR}/solver/umfpack_solver.cpp
  ${HERMES_COMMON_DIR}/solver/superlu.cpp
  ${HERMES_COMMON_DIR}/solver/precond_ml.cpp 
  ${HERMES_COMMON_DIR}/solver/precond_ifpack.cpp 
  ${HERMES_COMMON_DIR}/solver/eigensolver.cpp
  ${HERMES_COMMON_DIR}/solver/eigen.cpp
  ${HERMES_COMMON_DIR}/compat/fmemopen.cpp 
  ${HERMES_COMMON_DIR}/compat/c99_functions.cpp

  ${TEUCHOS_SOURCES}

  ${HERMES_COMMON_DIR}/python/python_api.cpp
  ${HERMES_COMMON_DIR}/python/python_engine.cpp
)

include(PickRealOrCplxLibs)
    
# Sets linker settings.
macro(LIBRARY_SETTINGS BIN)
#    target_link_libraries(${BIN} hermes_common)

  if(NOT MSVC)
      IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
          target_link_libraries(${BIN} "rt")
      ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      if(HAVE_TEUCHOS_BFD)
          # stacktrace code needs binutils:
          target_link_libraries(${BIN} bfd iberty)
      endif(HAVE_TEUCHOS_BFD)
  endif(NOT MSVC)
  
  # Link with all third party libraries corresponding to selected configuration options.
  # The respective variables will be empty for unselected TPL's.
  PICK_REAL_OR_CPLX_LIBS(${BIN} ${BIN})
  target_link_libraries(  ${BIN} 
      ${GLUT_LIBRARY} ${GLEW_LIBRARY}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${METIS_LIBRARY}
      ${UMFPACK_LIBRARIES}
      ${PARDISO_LIBRARY}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_LIBRARIES}
      ${SUPERLU_LIBRARY}
      ${LAPACK_LIBRARIES}
      ${PTHREAD_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
	  ${CLAPACK_LIBRARY} ${BLAS_LIBRARY} ${F2C_LIBRARY}
      ${ADDITIONAL_LIBS}        
      ${PYTHON_LIBRARY}
  )
endmacro(LIBRARY_SETTINGS)

# real version
if(H3D_REAL)
	add_library(${HERMES3D_REAL} SHARED ${SRC})
	set_target_properties(${HERMES3D_REAL} PROPERTIES COMPILE_FLAGS -DH3D_REAL)
	LIBRARY_SETTINGS(${HERMES3D_REAL})
endif(H3D_REAL)

# complex version
if(H3D_COMPLEX)
	add_library(${HERMES3D_CPLX} SHARED ${SRC})
	set_target_properties(${HERMES3D_CPLX} PROPERTIES COMPILE_FLAGS -DH3D_COMPLEX)
	LIBRARY_SETTINGS(${HERMES3D_CPLX})
endif(H3D_COMPLEX)

INCLUDE (CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES m)
CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)

# Create file with preprocessor definitions exposing the build settings to the source code.
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

#
# INSTALL
#

# install header files
file(GLOB INC           "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB INC_SPACE     "${CMAKE_CURRENT_SOURCE_DIR}/space/*.h")
file(GLOB INC_SHAPESET  "${CMAKE_CURRENT_SOURCE_DIR}/shapeset/*.h")
file(GLOB INC_LOADER    "${CMAKE_CURRENT_SOURCE_DIR}/loader/*.h")
file(GLOB INC_OUTPUT    "${CMAKE_CURRENT_SOURCE_DIR}/output/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")
file(GLOB INC_COMMON    "${HERMES_COMMON_DIR}/*.h")
file(GLOB INC_SOLVER    "${HERMES_COMMON_DIR}/solver/*.h")
file(GLOB INC_COMPAT    "${HERMES_COMMON_DIR}/compat/*.h")

install(FILES ${INC}           DESTINATION include/hermes3d/src)
install(FILES ${INC_SPACE}     DESTINATION include/hermes3d/src/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/hermes3d/src/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/hermes3d/src/loader)
install(FILES ${INC_OUTPUT}    DESTINATION include/hermes3d/src/output)
install(FILES ${INC_INTEGRALS} DESTINATION include/hermes3d/src/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/hermes3d/src/adapt)
install(FILES ${INC_COMMON}    DESTINATION include/hermes_common)
install(FILES ${INC_SOLVER}    DESTINATION include/hermes_common/solver)
install(FILES ${INC_COMPAT}    DESTINATION include/hermes_common/compat)

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
	set_target_properties(${LIB} PROPERTIES COMPILE_FLAGS ${FLAGS})
endmacro(BUILD_LIB_ALL)

# install library
if(H3D_REAL)
	if(MSVC)
		BUILD_LIB_ALL(${HERMES3D_REAL} "-D H3D_REAL -DEXPORT_HERMES_DLL")	
		
	else(MSVC)
		install(TARGETS ${HERMES3D_REAL} LIBRARY DESTINATION lib)
	endif(MSVC)
endif(H3D_REAL)

if(H3D_COMPLEX)
	if(MSVC)
		BUILD_LIB_ALL(${HERMES3D_CPLX} "-D H3D_COMPLEX -DEXPORT_HERMES_DLL")	
	else(MSVC)
		install(TARGETS ${HERMES3D_CPLX} LIBRARY DESTINATION lib)
	endif(MSVC)
endif(H3D_COMPLEX)
