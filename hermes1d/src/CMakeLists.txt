# Enabling multiprocessor build on MSVC
if(MSVC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif(MSVC)
set(SRC
    iterator.cpp qsort.cpp
    lobatto.cpp legendre.cpp
    discrete_problem.cpp solution.cpp space.cpp
    ogprojection.cpp
    ogprojection2.cpp
    linearizer.cpp quad_std.cpp transforms.cpp
    adapt.cpp graph.cpp h1_polys.cpp
    weakform.cpp
  ${HERMES_COMMON_DIR}/common.cpp
  ${HERMES_COMMON_DIR}/hermes_logging.cpp
  ${HERMES_COMMON_DIR}/callstack.cpp
  ${HERMES_COMMON_DIR}/error.cpp
  ${HERMES_COMMON_DIR}/utils.cpp
  ${HERMES_COMMON_DIR}/matrix.cpp
  ${HERMES_COMMON_DIR}/tables.cpp
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
  ${HERMES_COMMON_DIR}/compat/fmemopen.cpp 
  ${HERMES_COMMON_DIR}/compat/c99_functions.cpp
  ${HERMES_COMMON_DIR}/common_time_period.cpp
    )

add_definitions(-DCOMPLEX=std::complex<double>)

# Create file with preprocessor definitions exposing the build settings to the source code.
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)
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
  target_link_libraries(  ${BIN}
      ${GLUT_LIBRARY} ${GLEW_LIBRARY}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${METIS_LIBRARY}
      ${UMFPACK_LIBRARY}
      ${SUPERLU_LIBRARY}
      ${AMD_LIBRARY}
      ${PARDISO_LIBRARY}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_REAL_LIBRARIES}
      ${LAPACK_LIBRARIES}
      ${PTHREAD_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
	  ${CLAPACK_LIBRARY} ${BLAS_LIBRARY} ${F2C_LIBRARY}
      ${ADDITIONAL_LIBS}        
  )
endmacro(LIBRARY_SETTINGS)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

add_library(${HERMES1D} SHARED ${SRC})
if(DEBUG)
    set_target_properties(hermes1d-debug PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_DEBUG})
endif(DEBUG)
if(RELEASE)
    set_target_properties(hermes1d PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_RELEASE})
endif(RELEASE)
LIBRARY_SETTINGS(${HERMES1D})


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

#install
if(MSVC)
	BUILD_LIB_ALL(${HERMES1D} "-DH1D_REAL -DEXPORT_HERMES_DLL")
else(MSVC)
	install(TARGETS ${HERMES1D} RUNTIME DESTINATION bin ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
endif(MSVC)

file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes1d)

