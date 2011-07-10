project(hermes_common)

	# Linker settings.
	include(ConfigureRPATH)
		ConfigureRPATH()
		
	# Check for stacktrace and math library.
	include(CMakeCommon.txt)

	# Additional building steps required by MSVC.
	include(BuildAndInstallScripts)

	include_directories(${CMAKE_CURRENT_SOURCE_DIR})
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/solvers)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/python_API)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/third_party_codes/trilinos-teuchos)

	# Teuchos RCP
	if (NOT WITH_TRILINOS)
		set( TEUCHOS_SOURCES
			src/third_party_codes/trilinos-teuchos/Teuchos_TypeNameTraits.cpp
			src/third_party_codes/trilinos-teuchos/Teuchos_dyn_cast.cpp
			src/third_party_codes/trilinos-teuchos/Teuchos_Ptr.cpp
			src/third_party_codes/trilinos-teuchos/Teuchos_RCPNode.cpp
			src/third_party_codes/trilinos-teuchos/Teuchos_TestForException.cpp
			src/third_party_codes/trilinos-teuchos/Teuchos_TypeNameTraits.cpp
		)
	endif (NOT WITH_TRILINOS)

	if (WITH_PYTHON)
		set ( PYTHON_SOURCES
				src/python_API/python_api.cpp
				src/python_API/python_engine.cpp
			)
	endif (WITH_PYTHON)

	#
	# Source files for the Hermes_common library.
	#
	set(SRC
		src/hermes_logging.cpp
		src/common_time_period.cpp
		src/callstack.cpp
		src/error.cpp
		src/matrix.cpp
		src/tables.cpp
		src/qsort.cpp
		src/c99_functions.cpp
		src/ord.cpp
		
		src/third_party_codes/trilinos-teuchos/Teuchos_stacktrace.cpp
		
		src/solvers/dpinterface.cpp
		src/solvers/linear_solver.cpp
		src/solvers/nonlinear_solver.cpp
		src/solvers/nox_solver.cpp
		src/solvers/epetra.cpp
		src/solvers/aztecoo_solver.cpp
		src/solvers/amesos_solver.cpp
		src/solvers/mumps_solver.cpp
		src/solvers/superlu_solver.cpp
		src/solvers/superlu_real.cpp
		src/solvers/superlu_cplx.cpp
		src/solvers/petsc_solver.cpp
		src/solvers/umfpack_solver.cpp
		src/solvers/precond_ml.cpp
		src/solvers/precond_ifpack.cpp
	 # src/solvers/eigensolver.cpp
	 # src/solvers/eigen.cpp
		
		${TEUCHOS_SOURCES}
		${PYTHON_SOURCES}
	)

	#
	# Subdirectories with CMakeLists.txt
	#
	if(HERMES_COMMON_WITH_TESTS)
			enable_testing()
			add_subdirectory(tests)
	endif(HERMES_COMMON_WITH_TESTS)

	if(WITH_PYTHON)
			set(CYTHON_INCLUDE_DIRECTORIES -I${CMAKE_CURRENT_SOURCE_DIR})
			add_subdirectory(hermes_common_python)
	endif(WITH_PYTHON)

	# Create file with preprocessor definitions exposing the build settings to the source code.
	configure_file(
		config.h.cmake
		include/config.h
	)

	#
	# Build and install the library.
	#
	include(BuildAndInstallScripts)

	# This is needed e.g. for eigenvalue problems.
	# FIXME: Find which projects need it and enclose their CMake scripts in
	# if(WITH_PYTHON) ... endif(WITH_PYTHON). The following should then be not needed,
	# as it has been already done in ../CMakeLists.txt
	############ WORKING ON REPLACING THIS BY A C++ EIGENPROBLEM SOLVER.
	#if(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)
	#    find_package(Python REQUIRED)
	#    find_package(NumPy REQUIRED)
	#    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
	#endif(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)

	# Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
	# and set its build properties (compile and link flags, installation directories).
	macro(BUILD_COMMON_LIB HERMES_COMMON_LIB BUILD_FLAGS)
		add_library(${HERMES_COMMON_LIB} SHARED ${SRC})
		set_target_properties(${HERMES_COMMON_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})

		# Link with all third party libraries corresponding to selected configuration
		# options. Linking of optional libraries for the different Hermeses ({1,2,3}D)
		# is controled in their respective CMakeLists.txt. Note that for unselected
		# TPL's, the respective variables will be empty.

		target_link_libraries(  ${HERMES_COMMON_LIB}
				${EXODUSII_LIBRARIES}
				${HDF5_LIBRARY}
				${METIS_LIBRARY}
				${UMFPACK_LIBRARIES}
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

		if(NOT MSVC)
				IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
						target_link_libraries(${HERMES_COMMON_LIB} "rt")
				ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
				if(HAVE_TEUCHOS_BFD)
						# stacktrace code needs binutils:
						target_link_libraries(${HERMES_COMMON_LIB} bfd iberty)
				endif(HAVE_TEUCHOS_BFD)
		endif(NOT MSVC)
	endmacro(BUILD_COMMON_LIB)

	#
	# Configure building and installation of all versions of hermes_common needed by
	# the selected Hermes configurations.
	#
	if(MSVC)
		BUILD_COMMON_LIB(${HERMES_COMMON_LIB} ${DEBUG_FLAGS})
		ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB}
												 ${HERMES_COMMON_LIB_DEBUG}
												 ${HERMES_COMMON_LIB_RELEASE})
		INSTALL_LIB(${HERMES_COMMON_LIB})
	else(MSVC)
		if(HERMES_COMMON_DEBUG)
			set(BUILD_FLAGS "${DEBUG_FLAGS}")
			BUILD_COMMON_LIB(${HERMES_COMMON_LIB_DEBUG} ${BUILD_FLAGS})
			INSTALL_LIB(${HERMES_COMMON_LIB_DEBUG})
		endif(HERMES_COMMON_DEBUG)
		
		if(HERMES_COMMON_RELEASE)
			set(BUILD_FLAGS "${RELEASE_FLAGS}")
			BUILD_COMMON_LIB(${HERMES_COMMON_LIB_RELEASE} ${BUILD_FLAGS})
			INSTALL_LIB(${HERMES_COMMON_LIB_RELEASE})
		endif(HERMES_COMMON_RELEASE)
	endif(MSVC)

	# Install header files.

	file(GLOB INC_COMMON    "include/*.h")
	file(GLOB INC_SOLVER    "include/solvers/*.h")
	file(GLOB INC_PYTHON    "include/python_API/*.h")
	file(GLOB INC_TP        "include/third_party_codes/trilinos-teuchos/*.h*")
	install(FILES ${INC_COMMON}    DESTINATION ${TARGET_ROOT}/include/hermes_common)
	install(FILES ${INC_SOLVER}    DESTINATION ${TARGET_ROOT}/include/hermes_common/solvers)
	install(FILES ${INC_PYTHON}    DESTINATION ${TARGET_ROOT}/include/hermes_common/python_API)
	install(FILES ${INC_TP}        DESTINATION ${TARGET_ROOT}/include/hermes_common/third_party_codes/trilinos-teuchos)
