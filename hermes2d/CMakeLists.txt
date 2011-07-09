project(hermes2d)
	# Linker settings.
	include(ConfigureRPATH)
	ConfigureRPATH()
		
	IF(H2D_DEBUG)
		set(DEBUG YES)
	ELSE(H2D_DEBUG)
		set(DEBUG NO)
	ENDIF(H2D_DEBUG)

	IF(H2D_RELEASE)
		set(RELEASE YES)
	ELSE(H2D_RELEASE)
		set(RELEASE NO)
	ENDIF(H2D_RELEASE)

	#
	# Version configuration.
	#
	if (NOT DEBUG AND NOT RELEASE)
		message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified. 
												 Re-run with DEBUG and/or RELEASE set to YES.")
	endif (NOT DEBUG AND NOT RELEASE)

	if(MSVC)
		set(HERMES2D "hermes2d")              # Base name of the library.
		set(HERMES2D_RELEASE ${HERMES2D})     # Name of the release version.
		set(HERMES2D_DEBUG ${HERMES2D}-debug) # Name of the debug version.
			
	else(MSVC)
		if(DEBUG)
			set(HERMES2D_DEBUG "hermes2d-debug")
			set(HERMES2D ${HERMES2D_DEBUG})
			# set(CMAKE_BUILD_TYPE Debug) 
		# This does not work with multi-version configurations,
		# because it is linked to the project (i.e. we would have
		# to have separate projects for Debug and Release)
		endif(DEBUG)

		if(RELEASE)
			set(HERMES2D_RELEASE "hermes2d")
			set(HERMES2D ${HERMES2D_RELEASE})
			# set(CMAKE_BUILD_TYPE Release)
		endif(RELEASE)
	endif(MSVC)

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
	# Subdirectories.
	#
	add_subdirectory(src)

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
