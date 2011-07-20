project(hermes2d)
	# Determine which version (debug x release) to build.
	if (NOT H2D_DEBUG AND NOT H2D_RELEASE)
		message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified. 
												 Re-run with H2D_DEBUG and/or H2D_RELEASE set to YES.")
	endif (NOT H2D_DEBUG AND NOT H2D_RELEASE)

	if(MSVC)
		set(HERMES2D "hermes2d")              # Base name of the library.
		set(HERMES2D_RELEASE ${HERMES2D})     # Name of the release version.
		set(HERMES2D_DEBUG ${HERMES2D}-debug) # Name of the debug version.
			
	else(MSVC)
		if(H2D_DEBUG)
			set(HERMES2D_DEBUG "hermes2d-debug")
			set(HERMES2D ${HERMES2D_DEBUG})
			# set(CMAKE_BUILD_TYPE Debug) 
		# This does not work with multi-version configurations,
		# because it is linked to the project (i.e. we would have
		# to have separate projects for Debug and Release)
		endif(H2D_DEBUG)

		if(H2D_RELEASE)
			set(HERMES2D_RELEASE "hermes2d")
			set(HERMES2D ${HERMES2D_RELEASE})
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
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/include/solvers)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes_common/include/python_API)
	include_directories(${CMAKE_HOME_DIRECTORY}/hermes2d)
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
		src/forms.cpp
		src/range.cpp
		src/asmlist.cpp
		src/newton_solver.cpp
		src/picard_solver.cpp

		src/adapt/adapt.cpp
		src/adapt/kelly_type_adapt.cpp

		src/boundary_conditions/essential_boundary_conditions.cpp

		src/function/transformable.cpp
		src/function/function.cpp
		src/function/mesh_function.cpp
		src/function/exact_solution.cpp
		src/function/solution.cpp 
		src/function/filter.cpp
		src/function/hermes_function.cpp

		src/mesh/refmap.cpp 
		src/mesh/curved.cpp
		src/mesh/refinement_type.cpp 
		src/mesh/element_to_refine.cpp
		src/mesh/exodusii.cpp 
		src/mesh/hash.cpp 
		src/mesh/h2d_reader.cpp
		src/mesh/mesh.cpp
		src/mesh/traverse.cpp
		src/mesh/mesh_data.cpp

		src/quadrature/limit_order.cpp
		src/quadrature/quad_std.cpp

		src/refinement_selectors/selector.cpp 
		src/refinement_selectors/order_permutator.cpp 
		src/refinement_selectors/optimum_selector.cpp 
		src/refinement_selectors/proj_based_selector.cpp 
		src/refinement_selectors/l2_proj_based_selector.cpp 
		src/refinement_selectors/h1_proj_based_selector.cpp 
		src/refinement_selectors/hcurl_proj_based_selector.cpp

		src/shapeset/shapeset.cpp 
		src/shapeset/shapeset_h1_ortho.cpp 
		src/shapeset/shapeset_h1_jacobi.cpp 
		src/shapeset/shapeset_h1_quad.cpp
		src/shapeset/shapeset_hc_legendre.cpp 
		src/shapeset/shapeset_hc_gradleg.cpp
		src/shapeset/shapeset_hd_legendre.cpp
		src/shapeset/shapeset_l2_legendre.cpp
		src/shapeset/precalc.cpp 

		src/space/space.cpp 
		src/space/space_h1.cpp 
		src/space/space_hcurl.cpp 
		src/space/space_l2.cpp
		src/space/space_hdiv.cpp

		src/views/base_view.cpp 
		src/views/mesh_view.cpp 
		src/views/order_view.cpp 
		src/views/scalar_view.cpp 
		src/views/stream_view.cpp 
		src/views/vector_base_view.cpp 
		src/views/vector_view.cpp 
		src/views/view.cpp 
		src/views/view_data.cpp 
		src/views/view_support.cpp
		src/views/linearizer.cpp
		src/views/orderizer.cpp
		src/views/vectorizer.cpp

		src/weakform/weakform.cpp 

		src/neighbor.cpp
		src/graph.cpp
		src/ogprojection.cpp
		src/hermes2d_common_defs.cpp  
		src/discrete_problem.cpp
		src/runge_kutta.cpp
		src/spline.cpp

		src/weakform_library/weakforms_elasticity.cpp
		src/weakform_library/weakforms_h1.cpp
		src/weakform_library/weakforms_hcurl.cpp
		src/weakform_library/weakforms_maxwell.cpp
		src/weakform_library/weakforms_neutronics.cpp
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
    if(H2D_DEBUG)
      set(FLAGS "${DEBUG_FLAGS}")
    elseif(H2D_DEBUG)
      set(FLAGS "${RELEASE_FLAGS}")
    endif(H2D_DEBUG)

		BUILD_2D_LIB(${HERMES2D} ${HERMES_COMMON_LIB} ${FLAGS})
		ADD_MSVC_BUILD_FLAGS(${HERMES2D} ${HERMES2D_DEBUG} ${HERMES2D_RELEASE})
		INSTALL_LIB(${HERMES2D})
	else(MSVC)
		if(H2D_DEBUG)
			BUILD_2D_LIB(${HERMES2D_DEBUG} ${HERMES_COMMON_LIB_DEBUG} ${DEBUG_FLAGS})
			INSTALL_LIB(${HERMES2D_DEBUG})
		endif(H2D_DEBUG)
		if(H2D_RELEASE)
			BUILD_2D_LIB(${HERMES2D_RELEASE} ${HERMES_COMMON_LIB_RELEASE} ${RELEASE_FLAGS})
			INSTALL_LIB(${HERMES2D_RELEASE})
		endif(H2D_RELEASE)
	endif(MSVC)

	# Install header files.

	file(GLOB INC            "include/*.h")
	file(GLOB INC_FORM       "include/form/*.h")
	file(GLOB INC_SPACE      "include/space/*.h")
	file(GLOB INC_MESH       "include/mesh/*.h")
	file(GLOB INC_BC         "include/boundary_conditions/*.h")
	file(GLOB INC_SHAPESET   "include/shapeset/*.h")
	file(GLOB INC_FUNCTION   "include/function/*.h")
	file(GLOB INC_INTEGRALS  "include/integrals/*.h")
	file(GLOB INC_LINEARIZER "include/linearizer/*.h")
	file(GLOB INC_ADAPT      "include/adapt/*.h")
	file(GLOB INC_QUAD       "include/quadrature/*.h")
	file(GLOB INC_REF_SEL    "include/refinement_selectors/*.h")
	file(GLOB INC_VIEWS      "include/views/*.h")
	file(GLOB INC_WEAKFORM   "include/weakform/*.h")
	file(GLOB INC_WEAKFORM_LIBRARY "include/weakform_library/*.h")

	install(FILES ${INC}           DESTINATION ${TARGET_ROOT}/include/hermes2d)
	install(FILES ${INC_SPACE}     DESTINATION ${TARGET_ROOT}/include/hermes2d/space)
	install(FILES ${INC_MESH}      DESTINATION ${TARGET_ROOT}/include/hermes2d/mesh)
	install(FILES ${INC_BC}        DESTINATION ${TARGET_ROOT}/include/hermes2d/boundary_conditions)
	install(FILES ${INC_SHAPESET}  DESTINATION ${TARGET_ROOT}/include/hermes2d/shapeset)
	install(FILES ${INC_FUNCTION}  DESTINATION ${TARGET_ROOT}/include/hermes2d/function)
	install(FILES ${INC_INTEGRALS} DESTINATION ${TARGET_ROOT}/include/hermes2d/integrals)
	install(FILES ${INC_LINEARIZER} DESTINATION ${TARGET_ROOT}/include/hermes2d/linearizer)
	install(FILES ${INC_ADAPT}     DESTINATION ${TARGET_ROOT}/include/hermes2d/adapt)
	install(FILES ${INC_QUAD}      DESTINATION ${TARGET_ROOT}/include/hermes2d/quadrature)
	install(FILES ${INC_REF_SEL}   DESTINATION ${TARGET_ROOT}/include/hermes2d/refinement_selectors)
	install(FILES ${INC_VIEWS}     DESTINATION ${TARGET_ROOT}/include/hermes2d/views)
	install(FILES ${INC_WEAKFORM}  DESTINATION ${TARGET_ROOT}/include/hermes2d/weakform)
	install(FILES ${INC_WEAKFORM_LIBRARY}  DESTINATION ${TARGET_ROOT}/include/hermes2d/weakform_library)

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
