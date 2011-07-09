project(hermes2d)
	# Determine which version (debug x release) to build.
	if (NOT H2D_DEBUG AND NOT H2D_RELEASE)
		message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified. 
												 Re-run with H2D_DEBUG and/or H2D_RELEASE set to YES.")
	endif (NOT H2D_DEBUG AND NOT H2D_RELEASE)

	if(MSVC)
		set(HERMES2D "hermes2d")              # Base name of the library.
		set(HERMES2D_H2D_RELEASE ${HERMES2D})     # Name of the release version.
		set(HERMES2D_H2D_DEBUG ${HERMES2D}-debug) # Name of the debug version.
			
	else(MSVC)
		if(H2D_DEBUG)
			set(HERMES2D_H2D_DEBUG "hermes2d-debug")
			set(HERMES2D ${HERMES2D_H2D_DEBUG})
			# set(CMAKE_BUILD_TYPE Debug) 
		# This does not work with multi-version configurations,
		# because it is linked to the project (i.e. we would have
		# to have separate projects for Debug and Release)
		endif(H2D_DEBUG)

		if(H2D_RELEASE)
			set(HERMES2D_H2D_RELEASE "hermes2d")
			set(HERMES2D ${HERMES2D_H2D_RELEASE})
			# set(CMAKE_BUILD_TYPE Release)
		endif(H2D_RELEASE)
	endif(MSVC)
	
	# Linker settings.
	include(ConfigureRPATH)
	ConfigureRPATH()

	#
	# Optional libraries and settings for H2D.
	#
	if(H2D_WITH_GLUT)
		if(NOT APPLE)   # See hermes/OSX.cmake for APPLE-specific search procedures.
			find_package(GLUT REQUIRED)
			find_package(GLEW REQUIRED)
		endif(NOT APPLE)
	else(H2D_WITH_GLUT)
		add_definitions(-DNOGLUT)
	endif(H2D_WITH_GLUT)

	if(H2D_WITH_VIEWER_GUI)
		find_package(ANTTWEAKBAR REQUIRED)
		include_directories(${ANTTWEAKBAR_INCLUDE_DIR})
	endif(H2D_WITH_VIEWER_GUI)
	
	# Mesh format.
	if(WITH_HDF5)
		find_package(HDF5 REQUIRED)
		include_directories(${HDF5_INCLUDE_DIR})
	endif(WITH_HDF5)

	if(WITH_EXODUSII)
		find_package(EXODUSII REQUIRED)
		include_directories(${EXODUSII_INCLUDE_DIR})
	endif(WITH_EXODUSII)
	
	#
	# Header files location for the Hermes2d library.
	#
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/include)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/solvers/include)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/python_API/include)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/hermes2d/include)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/form)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/space)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/mesh)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/boundary_conditions)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/shapeset)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/function)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/integrals)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/linearizer)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/adapt)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/quadrature)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/refinement_selectors)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/views)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/weakform)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d/include/weakform_library)

  #
	# Source files for the Hermes2d library.
	#
	set(SRC
		forms.cpp
		range.cpp
		asmlist.cpp
		newton_solver.cpp
		picard_solver.cpp

		adapt/adapt.cpp
		adapt/kelly_type_adapt.cpp

		boundary_conditions/essential_boundary_conditions.cpp

		function/transformable.cpp
		function/function.cpp
		function/mesh_function.cpp
		function/exact_solution.cpp
		function/solution.cpp 
		function/filter.cpp
		function/hermes_function.cpp

		mesh/refmap.cpp 
		mesh/curved.cpp
		mesh/refinement_type.cpp 
		mesh/element_to_refine.cpp
		mesh/exodusii.cpp 
		mesh/hash.cpp 
		mesh/h2d_reader.cpp
		mesh/mesh.cpp
		mesh/traverse.cpp
		mesh/mesh_data.cpp

		quadrature/limit_order.cpp
		quadrature/quad_std.cpp

		refinement_selectors/selector.cpp 
		refinement_selectors/order_permutator.cpp 
		refinement_selectors/optimum_selector.cpp 
		refinement_selectors/proj_based_selector.cpp 
		refinement_selectors/l2_proj_based_selector.cpp 
		refinement_selectors/h1_proj_based_selector.cpp 
		refinement_selectors/hcurl_proj_based_selector.cpp

		shapeset/shapeset.cpp 
		shapeset/shapeset_h1_ortho.cpp 
		shapeset/shapeset_h1_jacobi.cpp 
		shapeset/shapeset_h1_quad.cpp
		shapeset/shapeset_hc_legendre.cpp 
		shapeset/shapeset_hc_gradleg.cpp
		shapeset/shapeset_hd_legendre.cpp
		shapeset/shapeset_l2_legendre.cpp
		shapeset/precalc.cpp 

		space/space.cpp 
		space/space_h1.cpp 
		space/space_hcurl.cpp 
		space/space_l2.cpp
		space/space_hdiv.cpp

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
		views/linearizer.cpp
		views/orderizer.cpp
		views/vectorizer.cpp

		weakform/weakform.cpp 

		neighbor.cpp
		graph.cpp
		ogprojection.cpp
		hermes2d_common_defs.cpp  
		discrete_problem.cpp
		runge_kutta.cpp
		spline.cpp

		weakform_library/weakforms_elasticity.cpp
		weakform_library/weakforms_h1.cpp
		weakform_library/weakforms_hcurl.cpp
		weakform_library/weakforms_maxwell.cpp
		weakform_library/weakforms_neutronics.cpp
	)
		
	#
	# Build and install the library.
	#
	include(BuildAndInstallScripts)

	# Add a target for a specified version of the Hermes library,
	# and set its build properties (compile and link flags, installation directories). Due to
	# different optional libraries for H1D, H2D and H3D, this macro is dimension-specific.
	macro(BUILD_2D_LIB HERMES_LIB HERMES_COMMON_LIB BUILD_FLAGS)
		add_library(${HERMES_LIB} SHARED ${SRC})

		set(_FLAGS ${BUILD_FLAGS})

		if(NOT REPORT_WITH_LOGO)
			set(_FLAGS "${_FLAGS} -DH2D_NO_LOGO")
		endif(NOT REPORT_WITH_LOGO)
		if(NOT REPORT_TO_FILE)
			set(_FLAGS "${_FLAGS} -DH2D_REPORT_NO_FILE")
		endif(NOT REPORT_TO_FILE)
		if(REPORT_RUNTIME_CONTROL)
			set(_FLAGS "${_FLAGS} -DH2D_REPORT_RUNTIME_CONTROL")
		endif(REPORT_RUNTIME_CONTROL)
		if(REPORT_ALL)
			set(_FLAGS "${_FLAGS} -DHERMES_REPORT_ALL")
		else(REPORT_ALL)
			if(REPORT_WARNING)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_WARNING")
			endif(REPORT_WARNING)
			if(REPORT_INTR_WARNING)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_INTR_WARNING")
			endif(REPORT_INTR_WARNING)
			if(REPORT_INFO)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_INFO")
			endif(REPORT_INFO)
			if(REPORT_VERBOSE)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_VERBOSE")
			endif(REPORT_VERBOSE)
			if(REPORT_TRACE)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_TRACE")
			endif(REPORT_TRACE)
			if(REPORT_TIME)
				set(_FLAGS "${_FLAGS} -DHERMES_REPORT_TIME")
			endif(REPORT_TIME)
		endif(REPORT_ALL)

		if(WITH_VIEWER_GUI)
			set(_FLAGS "${_FLAGS} -DENABLE_VIEWER_GUI")
		endif(WITH_VIEWER_GUI)

		set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${_FLAGS})

		target_link_libraries(  ${HERMES_LIB} 
			${HERMES_COMMON_LIB}
			${GLUT_LIBRARY} ${GLEW_LIBRARY}
			${ANTTWEAKBAR_LIBRARY}
			${LAPACK_LIBRARY}
			${CLAPACK_LIBRARY} ${BLAS_LIBRARY}
			${PYTHON_LIBRARY}
		)
	endmacro(BUILD_2D_LIB)

	# Build and install:
	if(MSVC)
		set(FLAGS "${MSVC_DEFINES}")
		BUILD_2D_LIB(${HERMES2D} ${HERMES_COMMON_LIB} ${FLAGS})
		ADD_MSVC_BUILD_FLAGS(${HERMES2D} ${HERMES2D_DEBUG} ${HERMES2D_RELEASE})
		INSTALL_LIB(${HERMES2D})
	else(MSVC)
	if(DEBUG)
		set(BUILD_FLAGS "${DEBUG_FLAGS}")
		BUILD_2D_LIB(${HERMES2D_DEBUG} ${HERMES_COMMON_LIB_DEBUG} ${BUILD_FLAGS})
		INSTALL_LIB(${HERMES2D_DEBUG})
	endif(DEBUG)
	if(RELEASE)
		set(BUILD_FLAGS "${RELEASE_FLAGS}")
		BUILD_2D_LIB(${HERMES2D_RELEASE} ${HERMES_COMMON_LIB_RELEASE} ${BUILD_FLAGS})
		INSTALL_LIB(${HERMES2D_RELEASE})
	endif(RELEASE)
	endif(MSVC)

	# Install header files.

	file(GLOB INC            "hermes2d/include/*.h")
	file(GLOB INC_FORM       "hermes2d/include/form/*.h")
	file(GLOB INC_SPACE      "hermes2d/include/space/*.h")
	file(GLOB INC_MESH       "hermes2d/include/mesh/*.h")
	file(GLOB INC_BC         "hermes2d/include/boundary_conditions/*.h")
	file(GLOB INC_SHAPESET   "hermes2d/include/shapeset/*.h")
	file(GLOB INC_FUNCTION   "hermes2d/include/function/*.h")
	file(GLOB INC_INTEGRALS  "hermes2d/include/integrals/*.h")
	file(GLOB INC_LINEARIZER "hermes2d/include/linearizer/*.h")
	file(GLOB INC_ADAPT      "hermes2d/include/adapt/*.h")
	file(GLOB INC_QUAD       "hermes2d/include/quadrature/*.h")
	file(GLOB INC_REF_SEL    "hermes2d/include/refinement_selectors/*.h")
	file(GLOB INC_VIEWS      "hermes2d/include/views/*.h")
	file(GLOB INC_WEAKFORM   "hermes2d/include/weakform/*.h")
	file(GLOB INC_WEAKFORM_LIBRARY "hermes2d/include/weakform_library/*.h")

	install(FILES ${INC}           DESTINATION include)
	install(FILES ${INC_FORM}      DESTINATION include/form)
	install(FILES ${INC_SPACE}     DESTINATION include/space)
	install(FILES ${INC_MESH}      DESTINATION include/mesh)
	install(FILES ${INC_BC}        DESTINATION include/boundaryconditions)
	install(FILES ${INC_SHAPESET}  DESTINATION include/shapeset)
	install(FILES ${INC_FUNCTION}  DESTINATION include/function)
	install(FILES ${INC_INTEGRALS} DESTINATION include/integrals)
	install(FILES ${INC_LINEARIZER} DESTINATION include/linearizer)
	install(FILES ${INC_ADAPT}     DESTINATION include/adapt)
	install(FILES ${INC_QUAD}      DESTINATION include/quadrature)
	install(FILES ${INC_REF_SEL}   DESTINATION include/refinement_selectors)
	install(FILES ${INC_VIEWS}     DESTINATION include/views)
	install(FILES ${INC_WEAKFORM}  DESTINATION include/weakform)
	install(FILES ${INC_WEAKFORM_LIBRARY}  DESTINATION include/weakform_library)

	if(H2D_WITH_TEST_EXAMPLES)
		enable_testing()
		add_subdirectory(test_examples)
	endif(H2D_WITH_TEST_EXAMPLES)

	if(H2D_WITH_TESTS)
		enable_testing()
		add_subdirectory(tests)
	endif(H2D_WITH_TESTS)
	
	#
	# Documentation
	#

	# This doesn't work yet:
	#add_subdirectory(doc)
	#add_custom_target(doc)
	#add_custom_command(
	#	SOURCE	${DOXYGEN_CONFIG_FILE}
	#	COMMAND	${DOXYGEN_BINARY}
	#	ARGS	${DOXYGEN_CONFIG_FILE}
	#	TARGET	doc
	#	OUTPUTS	${PROJECT_BINARY_DIR}/doc/html
	#)

	#add_custom_command(
	#	SOURCE	doc
	#	TARGET	doc
	#	DEPENDS	${PROJECT_BINARY_DIR}/doc/html
	#)
	message("\n")
