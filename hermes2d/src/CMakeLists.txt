# Enabling multiprocessor build on MSVC
if(MSVC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif(MSVC)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# source files for the hermes2d library
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
       neighbor.cpp
       numerical_flux.cpp
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
       h2d_common.cpp  
       weakform.cpp 
       discrete_problem.cpp
       forms.cpp
       mesh_parser.cpp mesh_lexer.cpp
       exodusii.cpp 
       h2d_reader.cpp
       runge_kutta.cpp
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
       ${HERMES_COMMON_DIR}/common.cpp       
       ${HERMES_COMMON_DIR}/logging.cpp    
       ${HERMES_COMMON_DIR}/hermes_logging.cpp
       ${HERMES_COMMON_DIR}/common_time_period.cpp
       ${HERMES_COMMON_DIR}/callstack.cpp
       ${HERMES_COMMON_DIR}/error.cpp
       ${HERMES_COMMON_DIR}/utils.cpp
       ${HERMES_COMMON_DIR}/matrix.cpp
       ${HERMES_COMMON_DIR}/tables.cpp
       ${HERMES_COMMON_DIR}/Teuchos_stacktrace.cpp
       ${HERMES_COMMON_DIR}/solver/nox.cpp 
       ${HERMES_COMMON_DIR}/solver/epetra.cpp 
       ${HERMES_COMMON_DIR}/solver/aztecoo.cpp 
       ${HERMES_COMMON_DIR}/solver/amesos.cpp
       ${HERMES_COMMON_DIR}/solver/mumps.cpp 
       ${HERMES_COMMON_DIR}/solver/pardiso.cpp 
       ${HERMES_COMMON_DIR}/solver/superlu.cpp
       ${HERMES_COMMON_DIR}/solver/petsc.cpp 
       ${HERMES_COMMON_DIR}/solver/umfpack_solver.cpp
       ${HERMES_COMMON_DIR}/solver/precond_ml.cpp 
       ${HERMES_COMMON_DIR}/solver/precond_ifpack.cpp 
       ${HERMES_COMMON_DIR}/compat/fmemopen.cpp 
       ${HERMES_COMMON_DIR}/compat/c99_functions.cpp
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
      ${ANTTWEAKBAR_LIBRARY}
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
      ${PYTHON_LIBRARIES}
      ${EXODUSII_LIBRARY}
  )
endmacro(LIBRARY_SETTINGS)

macro(BUILD_REPORT_FLAGS FLAGS)
set(REPORT_FLAGS "")
if(NOT REPORT_WITH_LOGO)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_NO_LOGO")
endif(NOT REPORT_WITH_LOGO)
if(NOT REPORT_TO_FILE)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_NO_FILE")
endif(NOT REPORT_TO_FILE)
if(REPORT_RUNTIME_CONTROL)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_RUNTIME_CONTROL")
endif(REPORT_RUNTIME_CONTROL)
if(REPORT_ALL)
	set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_ALL")
else(REPORT_ALL)
	if(REPORT_WARNING)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_WARNING")
	endif(REPORT_WARNING)
	if(REPORT_INTR_WARNING)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_INTR_WARNING")
	endif(REPORT_INTR_WARNING)
	if(REPORT_INFO)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_INFO")
	endif(REPORT_INFO)
	if(REPORT_VERBOSE)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_VERBOSE")
	endif(REPORT_VERBOSE)
	if(REPORT_TRACE)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_TRACE")
	endif(REPORT_TRACE)
	if(REPORT_TIME)
		set(REPORT_FLAGS "${REPORT_FLAGS} -DHERMES_REPORT_TIME")
	endif(REPORT_TIME)
endif(REPORT_ALL)
set(${FLAGS} "${REPORT_FLAGS}")
endmacro(BUILD_REPORT_FLAGS)

# Makes Win32 path from Unix-style patch which is used by CMAKE. Used when a path is provided to an OS utility.
macro(MAKE_PATH PATH_OUT PATH_IN)
	if(WIN32)
		string(REPLACE "/" "\\" ${PATH_OUT} ${PATH_IN})
	else(WIN32)
		set(${PATH_OUT} ${PATH_IN})
	endif(WIN32)
endmacro(MAKE_PATH)

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
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")
file(GLOB INC_REF_SEL   "${CMAKE_CURRENT_SOURCE_DIR}/ref_selectors/*.h")
file(GLOB INC_VIEWS     "${CMAKE_CURRENT_SOURCE_DIR}/views/*.h")
file(GLOB INC_COMMON    "${HERMES_COMMON_DIR}/*.h")
file(GLOB INC_SOLVER    "${HERMES_COMMON_DIR}/solver/*.h")
file(GLOB INC_COMPAT    "${HERMES_COMMON_DIR}/compat/*.h")

install(FILES ${INC}           DESTINATION include/hermes2d/src/)
install(FILES ${INC_SPACE}     DESTINATION include/hermes2d/src/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/hermes2d/src/shapeset)
install(FILES ${INC_INTEGRALS} DESTINATION include/hermes2d/src/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/hermes2d/src/adapt)
install(FILES ${INC_REF_SEL}   DESTINATION include/hermes2d/src/ref_selectors)
install(FILES ${INC_VIEWS}     DESTINATION include/hermes2d/src/views)
install(FILES ${INC_COMMON}    DESTINATION include/hermes_common)
install(FILES ${INC_SOLVER}    DESTINATION include/hermes_common/solver)
install(FILES ${INC_COMPAT}    DESTINATION include/hermes_common/compat)

# build and install MSVC version:
if(MSVC)
	set(MSVC_DLL_DEFINES "-DEXPORT_HERMES_DLL")
	if(H2D_REAL)
		BUILD_LIB_ALL(${HERMES_REAL_BIN} "${MSVC_DLL_DEFINES} -DH2D_REAL")
	endif(H2D_REAL)
	if(H2D_COMPLEX)
		BUILD_LIB_ALL(${HERMES_CPLX_BIN} "${MSVC_DLL_DEFINES} -DH2D_COMPLEX")
	endif(H2D_COMPLEX)

else(MSVC)
	set(CXX_OTHER_FLAGS "-Wuninitialized")

	# build and install real version:
	if(H2D_REAL)
		if(DEBUG)
			BUILD_LIB(hermes2d-real-debug "${CMAKE_CXX_FLAGS_DEBUG} -DH2D_REAL ${CXX_OTHER_FLAGS}")
		endif(DEBUG)
		if(RELEASE)
			BUILD_LIB(hermes2d-real "${CMAKE_CXX_FLAGS_RELEASE} -DH2D_REAL ${CXX_OTHER_FLAGS}")
		endif(RELEASE)
	endif(H2D_REAL)

	# build and install complex version:
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
