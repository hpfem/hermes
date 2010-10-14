# Enabling multiprocessor build on MSVC
if(MSVC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif(MSVC)


set(SRC
       hash.cpp 
       mesh.cpp 
       regul.cpp 
       refmap.cpp 
       curved.cpp
       transform.cpp 
       traverse.cpp
       limit_order.cpp
       precalc.cpp 
       solution.cpp 
       filter.cpp
       common_time_period.cpp
       ../common/callstack.cpp
       ../common/error.cpp
       ../common/timer.cpp
       ../common/utils.cpp
       ../common/Teuchos_stacktrace.cpp
       space/space.cpp 
       space/space_h1.cpp 
       space/space_hcurl.cpp 
       space/space_l2.cpp
       space/space_hdiv.cpp
       linear1.cpp 
       linear2.cpp 
       linear3.cpp 
       graph.cpp
       quad_std.cpp
       shapeset/shapeset.cpp 
       shapeset/shapeset_h1_ortho.cpp 
       shapeset/shapeset_h1_jacobi.cpp 
       shapeset/shapeset_h1_quad.cpp
       shapeset/shapeset_hc_legendre.cpp 
       shapeset/shapeset_hc_gradleg.cpp
       shapeset/shapeset_hd_legendre.cpp
       shapeset/shapeset_l2_legendre.cpp
       qsort.cpp norm.cpp
       trans.cpp
       ogprojection.cpp
       adapt/adapt.cpp
       refinement_type.cpp 
       element_to_refine.cpp
       ref_selectors/selector.cpp 
       ref_selectors/order_permutator.cpp 
       ref_selectors/optimum_selector.cpp 
       ref_selectors/proj_based_selector.cpp 
       ref_selectors/l2_proj_based_selector.cpp 
       ref_selectors/h1_proj_based_selector.cpp 
       ref_selectors/hcurl_proj_based_selector.cpp
       common.cpp 
       matrix.cpp 
       weakform.cpp 
       feproblem.cpp 
       solver/nox.cpp 
       solver/epetra.cpp 
       solver/aztecoo.cpp 
       solver/amesos.cpp
       solver/mumps.cpp 
       solver/pardiso.cpp 
       solver/petsc.cpp 
       solver/umfpack_solver.cpp
       solver/precond_ml.cpp 
       solver/precond_ifpack.cpp 
       forms.cpp
       mesh_parser.cpp mesh_lexer.cpp
       exodusii.cpp 
       h2d_reader.cpp
       views/base_view.cpp 
       views/mesh_view.cpp 
       views/order_view.cpp 
       views/scalar_view.cpp 
       views/stream_view.cpp 
       views/vector_base_view.cpp 
       views/vector_view.cpp 
       views/view.cpp 
       views/view_data.cpp 
       views/view_support.cpp
       compat/fmemopen.cpp 
       compat/c99_functions.cpp
       )

# Makes Win32 path from Unix-style patch which is used by CMAKE. Used when a path is provided to an OS utility.
macro(MAKE_PATH PATH_OUT PATH_IN)
	if(WIN32)
		string(REPLACE "/" "\\" ${PATH_OUT} ${PATH_IN})
	else(WIN32)
		set(${PATH_OUT} ${PATH_IN})
	endif(WIN32)
endmacro(MAKE_PATH)

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
      ${JUDY_LIBRARY}
      ${GLUT_LIBRARY} ${GLEW_LIBRARY}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${ANTTWEAKBAR_LIBRARY}
      ${UMFPACK_LIBRARIES}
      ${PARDISO_LIBRARY}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      ${PTHREAD_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
      ${ADDITIONAL_LIBS}        
  )
endmacro(LIBRARY_SETTINGS)

macro(BUILD_REPORT_FLAGS FLAGS)
set(REPORT_FLAGS "")
if(NOT REPORT_WITH_LOGO)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_NO_LOGO")
endif(NOT REPORT_WITH_LOGO)
if(NOT REPORT_TO_FILE)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_NO_FILE")
endif(NOT REPORT_TO_FILE)
if(REPORT_RUNTIME_CONTROL)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_RUNTIME_CONTROL")
endif(REPORT_RUNTIME_CONTROL)
if(REPORT_ALL)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_ALL")
else(REPORT_ALL)
	if(REPORT_WARNING)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_WARNING")
	endif(REPORT_WARNING)
	if(REPORT_INTR_WARNING)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_INTR_WARNING")
	endif(REPORT_INTR_WARNING)
	if(REPORT_INFO)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_INFO")
	endif(REPORT_INFO)
	if(REPORT_VERBOSE)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_VERBOSE")
	endif(REPORT_VERBOSE)
	if(REPORT_TRACE)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_TRACE")
	endif(REPORT_TRACE)
	if(REPORT_TIME)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DH2D_REPORT_TIME")
	endif(REPORT_TIME)
endif(REPORT_ALL)
set(${FLAGS} "${REPORT_FLAGS}")
endmacro(BUILD_REPORT_FLAGS)

macro(BUILD_LIB LIB FLAGS)
  BUILD_REPORT_FLAGS(REPORT_FLAGS)
	set(BUILD_FLAGS "${FLAGS} ${REPORT_FLAGS}")
	if(WITH_VIEWER_GUI)
		set(BUILD_FLAGS "${BUILD_FLAGS} -DENABLE_VIEWER_GUI")
	endif(WITH_VIEWER_GUI)

	add_library(${LIB} SHARED ${SRC})
	set_target_properties(${LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
	LIBRARY_SETTINGS(${LIB})
	install(TARGETS ${LIB}
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION lib
	)
endmacro(BUILD_LIB)

macro(BUILD_LIB_ALL LIB FLAGS)
	BUILD_REPORT_FLAGS(REPORT_FLAGS)
	set(BUILD_FLAGS "${FLAGS} ${REPORT_FLAGS}")
	if(WITH_VIEWER_GUI)
		set(BUILD_FLAGS "${BUILD_FLAGS} -DENABLE_VIEWER_GUI")
	endif(WITH_VIEWER_GUI)

	add_library(${LIB} SHARED ${SRC})
	set_target_properties(${LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
	set_target_properties(${LIB} PROPERTIES DEBUG_OUTPUT_NAME "${LIB}-debug")
	set_target_properties(${LIB} PROPERTIES RELEASE_OUTPUT_NAME ${LIB})
	LIBRARY_SETTINGS(${LIB})	
	
	install(TARGETS ${LIB}
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION lib)
		
	# Copy target to a place which is accessible by the system. This is crucial for successful debugging using MSVC.
	if(MSVC)
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
		unset(SOURCE_FILE)
	endif(MSVC)
endmacro(BUILD_LIB_ALL)

INCLUDE (CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES m)
CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)

# Create a file with preprocessor definitions (mainly to determine which TPL's 
# are present on the system during compilation):
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
                ${CMAKE_CURRENT_SOURCE_DIR}/config.h  )


# build MSVC version:
if(MSVC)
	set(MSVC_DLL_DEFINES "-D_HERMESDLL")
	if(H2D_REAL)
		BUILD_LIB_ALL(${HERMES_REAL_BIN} "${MSVC_DLL_DEFINES}")
	endif(H2D_REAL)
	if(H2D_COMPLEX)
		BUILD_LIB_ALL(${HERMES_CPLX_BIN} "${MSVC_DLL_DEFINES} -DH2D_COMPLEX")
	endif(H2D_COMPLEX)

else(MSVC)
	set(CXX_OTHER_FLAGS "-Wuninitialized")

	# build real version:
	if(H2D_REAL)
		if(DEBUG)
			BUILD_LIB(hermes2d-real-debug "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_OTHER_FLAGS}")
		endif(DEBUG)
		if(RELEASE)
			BUILD_LIB(hermes2d-real "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_OTHER_FLAGS}")
		endif(RELEASE)
	endif(H2D_REAL)

	# build complex version:
	if(H2D_COMPLEX)
		if(DEBUG)
			BUILD_LIB(hermes2d-cplx-debug "${CMAKE_CXX_FLAGS_DEBUG} -DH2D_COMPLEX ${CXX_OTHER_FLAGS}")
		endif(DEBUG)
		if(RELEASE)
			BUILD_LIB(hermes2d-cplx "${CMAKE_CXX_FLAGS_RELEASE} -DH2D_COMPLEX ${CXX_OTHER_FLAGS}")
		endif(RELEASE)
	endif(H2D_COMPLEX)

	if(DEBUG)
		add_custom_target(debug)
		add_dependencies(debug
			hermes2d-real-debug
			hermes2d-cplx-debug)
	endif(DEBUG)

	if(RELEASE)
		add_custom_target(release)
		add_dependencies(release
			hermes2d-real
			hermes2d-cplx)
	endif(RELEASE)
endif(MSVC)

file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d)
file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/compat/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d/compat)
file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/views/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d/views)
file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/ref_selectors/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d/ref_selectors)
