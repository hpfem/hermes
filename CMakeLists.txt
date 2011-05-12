set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/hermes_common/cmake/cxx_flag_overrides.cmake)
project(hermes)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/hermes_common)
set(HERMES_COMMON_LIB_REAL  hermes_common-real)
set(HERMES_COMMON_LIB_CPLX  hermes_common-cplx)


#ADDED BY COLMAN
IF (APPLE)
    IF (NOT CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
        IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            ## It's Mac OS X.

            SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Release is default." FORCE)

            SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING
                "compiling for OSX only in i386 format. Library with the
                 necessary symbol '_glutSetOption' was only available
                 in /usr/lib/libglut.dylib, which happened to be a pure
                 i386 library. Maybe fixed later." FORCE)
            IF (${CMAKE_VERBOSE_MAKEFILE} MATCHES "TRUE")
                MESSAGE(STATUS)
                MESSAGE(STATUS "*******************************************************************************")
                MESSAGE(STATUS "*        SUMMARY OF USED VARIABLES -> Details in AllVariables.txt             *")
                MESSAGE(STATUS "*******************************************************************************")
                MESSAGE(STATUS "COMPILING FOR OSX, i386 platform")
                MESSAGE(STATUS "CMAKE_ROOT ........................... = ${CMAKE_ROOT}")
                MESSAGE(STATUS "CMAKE_HOME_DIRECTORY ................. = ${CMAKE_HOME_DIRECTORY}")
                MESSAGE(STATUS "CMAKE_VERSION ........................ = ${CMAKE_VERSION}")
                MESSAGE(STATUS "CMAKE_GENERATOR ...................... = ${CMAKE_GENERATOR}")
                MESSAGE(STATUS "CMAKE_MODULE_PATH .................... = ${CMAKE_MODULE_PATH}")
                MESSAGE(STATUS "CMAKE_HOST_SYSTEM .................... = ${CMAKE_HOST_SYSTEM}")
                MESSAGE(STATUS "CMAKE_HOST_SYSTEM_VERSION ............ = ${CMAKE_HOST_SYSTEM_VERSION}")
                MESSAGE(STATUS "CMAKE_HOST_SYSTEM_NAME ............... = ${CMAKE_HOST_SYSTEM_NAME}")
                MESSAGE(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR .......... = ${CMAKE_HOST_SYSTEM_PROCESSOR}")
                MESSAGE(STATUS "CMAKE_HOST_UNIX ...................... = ${CMAKE_HOST_UNIX}")
                MESSAGE(STATUS "CMAKE_BUILD_TYPE ..................... = ${CMAKE_BUILD_TYPE}")
                MESSAGE(STATUS "CMAKE_C_COMPILER ..................... = ${CMAKE_C_COMPILER}")
                MESSAGE(STATUS "CMAKE_CXX_COMPILER ................... = ${CMAKE_CXX_COMPILER}")
            ENDIF()
        ENDIF ()
    ENDIF (NOT CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
ENDIF (APPLE)


#------------------------------------------------------------------------------
# Main configuration options.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### BEGIN DEFAULT CONFIGURATION ###

# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

# H1D:
set(WITH_H1D                  YES)
	# For use in the common code.
	set(H1D_REAL              YES)
	# Release and debug versions.
	set(H1D_DEBUG             YES)
	set(H1D_RELEASE           NO)
	# Optional parts of the library.
	set(H1D_WITH_EXAMPLES     YES)
	set(H1D_WITH_BENCHMARKS   YES)
	set(H1D_WITH_TESTS        YES)

# H2D:
set(WITH_H2D                  YES)
	# Real/complex version of the library (you can build either one or both):
	set(H2D_REAL              YES)
	set(H2D_COMPLEX           YES)
	# Release and debug versions.
	set(H2D_DEBUG             YES)
	set(H2D_RELEASE           YES)
	# Optional parts of the library.
	set(H2D_WITH_EXAMPLES     YES)
	set(H2D_WITH_BENCHMARKS   YES)
	set(H2D_WITH_TESTS        YES)
	set(H2D_WITH_TUTORIAL     YES)
	set(H2D_WITH_GLUT         YES)
	set(H2D_WITH_VIEWER_GUI   NO)
	set(H2D_WITH_UTIL         YES)

# H3D:		
set(WITH_H3D                  YES)
	# Real/complex version of the library.
	set(H3D_REAL              YES)
	set(H3D_COMPLEX           YES)
	# Release and debug versions.
	set(H3D_DEBUG             YES)
	set(H3D_RELEASE           NO)
	# Debugging options.
	set(H3D_DEBUG_ORDER       NO)
	set(H3D_PROFILING         NO)
	# Supported element types.
	set(WITH_TETRA            YES)
	set(WITH_HEX              YES)
	set(WITH_PRISM            NO)
	
	# Optional parts of the library.
	set(H3D_WITH_EXAMPLES     YES)
	set(H3D_WITH_BENCHMARKS   YES)
	set(H3D_WITH_TESTS        YES)
	set(H3D_WITH_TOOLS        NO)
	set(H3D_WITH_DEV_TESTS    NO)
	
# Hermes_common (linked to all H1D, H2D and H3D):
set(HERMES_COMMON_WITH_TESTS  YES)
set(HERMES_COMMON_CHECK_BOUNDARY_CONDITIONS  YES)

# Reporting and logging:
set(REPORT_WITH_LOGO        YES)  #logo will be shown
set(REPORT_TO_FILE          YES)  #reports will be mirrored to a file
set(REPORT_RUNTIME_CONTROL  NO)   #reporting will be controlled by preprocessor directives. 
                                  #The runtime control is intended just for Python wrappers.
set(REPORT_ALL              NO)   #do not report all events, just selected ones
set(REPORT_WARNING          YES)  #warnings will be reported
set(REPORT_INTR_WARNING     NO)   #warnings about itegration points will NOT be reported
set(REPORT_INFO             YES)  #info will be reported
set(REPORT_VERBOSE          NO)   #info details will not be reported
set(REPORT_TRACE            NO)   #code execution tracing will not be reported
set(REPORT_TIME             NO)   #time will not be measured and time measurement will not be reported
#set(REPORT_DEBUG           NO)   #debug events will depend on version which is compiled

#### Solvers ###
# Standard UMFPACK.
set(WITH_UMFPACK            YES)

# First option enables sequential SuperLU, which should reside in SUPERLU_ROOT/lib
# and SUPERLU_ROOT/include. If SUPERLU_ROOT is not provided, the library will be searched
# for in system directories. If SUPERLU_MT is set to YES, the multithreaded version of
# SuperLU will be used. In this case, there should be either superlu_mt_PTHREAD.a and/or 
# superlu_mt_OPENMP.a libraries present in the search path (the former will be chosen if 
# WITH_OPENMP is set to NO, the latter otherwise). There should also be a subdirectory
# SLU_MT in the include directory since SuperLU_MT uses different versions of the same-named
# header files as the serial library. 
#
set(WITH_SUPERLU            NO)
set(SUPERLU_ROOT            /opt/SuperLU_4.0/)
set(SUPERLU_MT              NO)

# If PETSc is enabled, PETSC_ROOT and PETSC_ARCH must be set as well;
# if PETSc is disabled, PETSC_ROOT and PETSC_ARCH will be ignored.
#
set(WITH_PETSC              NO)
set(PETSC_ROOT              /opt/petsc/petsc-3.1-p7)
set(PETSC_ARCH              petsc-arch-without-scalar-type)

# Sequential version of MUMPS may be enabled by setting YES on the following line
# and pointing to the package installation directory with MUMPS_ROOT. To enable the 
# parallel version, setting WITH_MPI to YES in addition to that should suffice. If link
# problems arise, it may be that the automatically found SCALAPACK libraries (reported
# during the CMake configuration phase) do not match those used when compiling MUMPS.
# In this case, set SCALAPACK_LIBRARIES manually in CMake.vars and reconfigure Hermes.
#
set(WITH_MUMPS              NO)
# set(MUMPS_ROOT            root/dir/of/mumps/)
# set(SCALAPACK_LIBRARIES   -lscalapack -lblacs)

# Enable support for Trilinos solvers.
set(WITH_TRILINOS           NO)


### Mesh format ###
set(WITH_EXODUSII           NO)
set(WITH_HDF5               NO)

### Others ###
# Parallel execution (tells the linker to use parallel versions of the selected 
# solvers, if available):
#
# Enable OpenMP for shared memory processing.
#  
set(WITH_OPENMP             NO)

# If MPI is enabled, the MPI library installed on the system should be found by 
# CMake automatically. If the found library doesn't match the one used to compile the
# particular MPI-dependent package, the other two options should be used to specify it.
#
set(WITH_MPI                NO)
# set(MPI_LIBRARIES         -lmpi)
# set(MPI_INCLUDE_PATH      /usr/include/openmpi

# Enable or disable Python.
set(WITH_PYTHON             YES)

# Set overrides:
if(WITH_PYTHON)
  set(REPORT_RUNTIME_CONTROL YES) #reuired by python wrappers
endif(WITH_PYTHON)

### end DEFAULT CONFIGURATION ###

# Allow to override the default values in CMake.vars:
include(CMake.vars OPTIONAL)

# Ensure that the correct version of hermes_common required by the selected
# Hermes libraries will be built. Note that the debug/release decision has to
# be made also here and not in ./hermes_common, so as to expose 
# HERMES_COMMON_LIB_{REAL/CPLX}_{DEBUG/RELEASE} to all Hermeses. The logic is
# as follows:
# If any HXD (X={1/2/3}) is configured for a Debug version, Debug version of 
# hermes_common will be built (the same for Release).  
# If any HXD is configured for a Complex version, Complex version of 
# hermes_common will be built (the same for Real).  

macro (DETERMINE_HERMES_COMMON_REAL_COMPLEX HXD)
    if (WITH_${HXD})
        if (${HXD}_REAL)
            set(HERMES_COMMON_REAL  YES)
        endif (${HXD}_REAL)
        if (${HXD}_COMPLEX)
            set(HERMES_COMMON_COMPLEX YES)
        endif (${HXD}_COMPLEX)
    endif (WITH_${HXD})   
endmacro(DETERMINE_HERMES_COMMON_REAL_COMPLEX)

set(HERMES_COMMON_REAL    NO)
set(HERMES_COMMON_COMPLEX NO)

DETERMINE_HERMES_COMMON_REAL_COMPLEX(H1D)
DETERMINE_HERMES_COMMON_REAL_COMPLEX(H2D)
DETERMINE_HERMES_COMMON_REAL_COMPLEX(H3D)

macro (DETERMINE_HERMES_COMMON_RELEASE_DEBUG HXD)
    if (WITH_${HXD})
        if(MSVC)
            # In MSVC, Release/Debug is distinguished in a special way
            # (see hermes_common/cmake/BuildAndInstallScripts.cmake)
            if (HERMES_COMMON_REAL)
                set(HERMES_COMMON_LIB_REAL_DEBUG ${HERMES_COMMON_LIB_REAL}-debug)
                set(HERMES_COMMON_LIB_REAL_RELEASE ${HERMES_COMMON_LIB_REAL})
            endif (HERMES_COMMON_REAL)
            if (HERMES_COMMON_COMPLEX)
                set(HERMES_COMMON_LIB_CPLX_RELEASE ${HERMES_COMMON_LIB_CPLX})
                set(HERMES_COMMON_LIB_CPLX_DEBUG ${HERMES_COMMON_LIB_CPLX}-debug)
            endif (HERMES_COMMON_COMPLEX)
        else(MSVC)        
          if(${HXD}_DEBUG)
            set(HERMES_COMMON_DEBUG YES)
            
            if(HERMES_COMMON_REAL)
              set(HERMES_COMMON_LIB_REAL_DEBUG ${HERMES_COMMON_LIB_REAL}-debug)
            endif(HERMES_COMMON_REAL)
            if(HERMES_COMMON_COMPLEX)
              set(HERMES_COMMON_LIB_CPLX_DEBUG  ${HERMES_COMMON_LIB_CPLX}-debug)
            endif(HERMES_COMMON_COMPLEX)
          endif(${HXD}_DEBUG)

          if(${HXD}_RELEASE)
            set(HERMES_COMMON_RELEASE YES)
            
            if(HERMES_COMMON_REAL)
              set(HERMES_COMMON_LIB_REAL_RELEASE ${HERMES_COMMON_LIB_REAL})
            endif(HERMES_COMMON_REAL)
            if(HERMES_COMMON_COMPLEX)
              set(HERMES_COMMON_LIB_CPLX_RELEASE ${HERMES_COMMON_LIB_CPLX})
            endif(HERMES_COMMON_COMPLEX)
          endif(${HXD}_RELEASE)
        endif(MSVC)
    endif (WITH_${HXD})   
endmacro(DETERMINE_HERMES_COMMON_RELEASE_DEBUG)

if (NOT MSVC)
  set(HERMES_COMMON_RELEASE   NO)
  set(HERMES_COMMON_DEBUG     NO)
endif (NOT MSVC)

DETERMINE_HERMES_COMMON_RELEASE_DEBUG(H1D)
DETERMINE_HERMES_COMMON_RELEASE_DEBUG(H2D)
DETERMINE_HERMES_COMMON_RELEASE_DEBUG(H3D)

#------------------------------------------------------------------------------
# Find necessary packages according to the settings above.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_package(PTHREAD REQUIRED)

# Python
if(WITH_PYTHON)
  find_package(Python REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(WITH_PYTHON)

# Solvers (TODO: This could maybe go into hermes_common/solver/CMakeLists.txt).
if(MSVC)
	if(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
	    # The search for CLAPACK and associated libraries will be performed 
	    # every time the project is configured in MSVC (TODO: Find out if 
	    # CLAPACK_LIBRARY, BLAS_LIBRARY, F2C_LIBRARY and CLAPACK_INCLUDE_DIRS
	    # could be cached to improve subsequent configuration performance).
	    find_package(CLAPACK REQUIRED)
	    include_directories(${CLAPACK_INCLUDE_DIRS})
	endif(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
else(MSVC)
  # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
  # FindLAPACK is Fortran-only so Fortran must be enabled as well.
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    if (NOT LAPACK_FOUND)
        enable_language(Fortran)
        find_package(LAPACK REQUIRED)
            
        # If no error occured, LAPACK library has been found. Save the path to
        # it to cache, so that it will not be searched for during next 'cmake .'
        set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES} 
                CACHE STRING  "Path to LAPACK/BLAS libraries.")
        set(LAPACK_FOUND      YES 
                CACHE STRING  "Have LAPACK/BLAS libraries been found?")
    endif (NOT LAPACK_FOUND)
    
    add_definitions(-DWITH_BLAS)	
  endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
endif(MSVC)  

if(WITH_PETSC)
	find_package(PETSC REQUIRED)
	# For some header files, PETSc uses the same name for both real and complex 
	# versions (they only reside in different directories), so we can't globally
	# set the architecture-specific include directory here, but rather on a 
	# per-target basis.
	include_directories(${COMMON_PETSC_INCLUDE_DIRS})
endif(WITH_PETSC)

if(WITH_UMFPACK)
	find_package(UMFPACK REQUIRED)
	include_directories(${UMFPACK_INCLUDE_DIRS})
endif(WITH_UMFPACK)

if(WITH_MUMPS)
	find_package(MUMPS REQUIRED)
	if(WITH_MPI)
		find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS 
  endif(WITH_MPI)
  include_directories(${MUMPS_INCLUDE_DIR})
endif(WITH_MUMPS)

if(WITH_TRILINOS)
	find_package(TRILINOS REQUIRED)
	include_directories(${TRILINOS_INCLUDE_DIR})
endif(WITH_TRILINOS)

if(WITH_SUPERLU)
  find_package(SUPERLU REQUIRED)
  include_directories(${SUPERLU_INCLUDE_DIR})
  if(SUPERLU_MT)
    add_definitions(-DSLU_MT)
  endif(SUPERLU_MT)
endif(WITH_SUPERLU)

# SMP.
if(WITH_OPENMP)
	find_package(OpenMP REQUIRED)
	set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif(WITH_OPENMP)

# Mesh format.
if(WITH_HDF5)
	find_package(HDF5 REQUIRED)
	include_directories(${HDF5_INCLUDE_DIR})
endif(WITH_HDF5)

if(WITH_EXODUSII)
	find_package(EXODUSII REQUIRED)
	include_directories(${EXODUSII_INCLUDE_DIR})
endif(WITH_EXODUSII)

# If using any package that requires MPI (e.g. parallel versions of MUMPS, PETSC).
if(WITH_MPI)
  if(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH) # If MPI was not defined by the user
    find_package(MPI REQUIRED)
  else(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
    if(NOT EXISTS MPI_INCLUDE_PATH)
      message(FATAL_ERROR "MPI include directory ${MPI_INCLUDE_PATH} not found")
    endif(NOT EXISTS MPI_INCLUDE_PATH)
  endif(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
  include_directories(MPI_INCLUDE_PATH)	  
endif(WITH_MPI)

# Use modified search procedures for some libraries on MacOS.
IF (APPLE)
    include(${CMAKE_HOME_DIRECTORY}/OSX.cmake)
ENDIF (APPLE)


#------------------------------------------------------------------------------
# Project-wide compiler settings.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Enable debugging symbols for all files all the time:
# set(CMAKE_CXX_FLAGS "-g ${CMAKE_CXX_FLAGS}")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Stop on first error (sometimes gcc doesn't print all the possible
    # candidates when a function signature is wrong if this is enabled, so we
    # leave it disabled by default):
    # set(CMAKE_CXX_FLAGS "-Wfatal-errors ${CMAKE_CXX_FLAGS}")

    # All warnings are errors
    #set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")

    # Disable the "unused variable/parameters/functions" warnings (as those are
    # not really problematic, it doesn't cause any harm):
    set(CMAKE_CXX_FLAGS "-Wno-unused-parameter -Wno-unused-variable -Wno-unused-function ${CMAKE_CXX_FLAGS}")

    # Disable all switch warnings for now (until all code is fixed)
    set(CMAKE_CXX_FLAGS "-Wno-switch-enum -Wno-switch -Wno-switch-default ${CMAKE_CXX_FLAGS}")

    # Enable all warnings:
    set(CMAKE_CXX_FLAGS "-W -Wuninitialized -Wall ${CMAKE_CXX_FLAGS}")

    set(RELEASE_FLAGS "-DNDEBUG -O3")
    set(DEBUG_FLAGS   "-g")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

# Enabling multiprocessor build on MSVC
if(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") 
endif(MSVC)


#------------------------------------------------------------------------------
# Project-wide linker settings.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include(ConfigureRPATH)
ConfigureRPATH()

#------------------------------------------------------------------------------
# Add project subdirectories and custom targets.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

add_subdirectory(${HERMES_COMMON_DIR})
include_directories(${HERMES_COMMON_DIR})

if(WITH_H1D)
    add_subdirectory(hermes1d)
endif(WITH_H1D)

if(WITH_H2D)
    add_subdirectory(hermes2d)
endif(WITH_H2D)

if(WITH_H3D)
    add_subdirectory(hermes3d)
endif(WITH_H3D)

if(H1D_WITH_TESTS OR H2D_WITH_TESTS OR H3D_WITH_TESTS OR HERMES_COMMON_WITH_TESTS)
  enable_testing()
  add_custom_target(test COMMAND ctest)
  add_custom_target(test-quick COMMAND ctest -LE slow -j9)
endif(H1D_WITH_TESTS OR H2D_WITH_TESTS OR H3D_WITH_TESTS OR HERMES_COMMON_WITH_TESTS)


#------------------------------------------------------------------------------
# Report.
#~~~~~~~~

message("\nConfiguration results")
message("---------------------")
message("Build with Python: ${WITH_PYTHON}")
message("Build with UMFPACK: ${WITH_UMFPACK}")
message("Build with PETSC: ${WITH_PETSC}")
message("Build with MUMPS: ${WITH_MUMPS}")
message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
message("Build with TRILINOS: ${WITH_TRILINOS}")
message("Build with MPI: ${WITH_MPI}")
message("Build with OPENMP: ${WITH_OPENMP}")
if(HAVE_TEUCHOS_STACKTRACE)
    message("Print Teuchos stacktrace on segfault: YES")
else(HAVE_TEUCHOS_STACKTRACE)
    message("Print Teuchos stacktrace on segfault: NO")
endif(HAVE_TEUCHOS_STACKTRACE)
if(REPORT_RUNTIME_CONTROL)
  if(WITH_PYTHON)
    message("Report controlled by: RUNTIME and PREP. DIRECTIVES (enforced by Python)")
  else(WITH_PYTHON)
    message("Report controlled by: RUNTIME and PREP. DIRECTIVES (enforced by Python)")
  endif(WITH_PYTHON)
else(REPORT_RUNTIME_CONTROL)
  message("Report controlled by: PREPROCESSOR DIRECTIVES")
endif(REPORT_RUNTIME_CONTROL)
message("Report with logo: ${REPORT_WITH_LOGO}")
if(REPORT_ALL)
  message("Report all events: YES")
else(REPORT_ALL)
  message("Report warnings: ${REPORT_WARNING}")
  message("Report integration warnings: ${REPORT_INTR_WARNING}")
  message("Report results/info: ${REPORT_INFO}")
  message("Report results/info details: ${REPORT_VERBOSE}")
  message("Report program execution: ${REPORT_TRACE}")
  message("Report internal time measurements: ${REPORT_TIME}")
endif(REPORT_ALL)
message("---------------------")
message("Hermes common library:")
message("\tBuild Real version: ${HERMES_COMMON_REAL}")
message("\tBuild Complex version: ${HERMES_COMMON_COMPLEX}")
message("\tBuild Debug version: ${HERMES_COMMON_DEBUG}")
message("\tBuild Release version: ${HERMES_COMMON_RELEASE}")
message("\tBuild with tests: ${HERMES_COMMON_WITH_TESTS}")
message("\tBuild with check boundary conditions: ${HERMES_COMMON_CHECK_BOUNDARY_CONDITIONS}")

message("Build Hermes1D: ${WITH_H1D}")
if(WITH_H1D)
	message("\tBuild Hermes1D Debug version: ${H1D_DEBUG}")
	message("\tBuild Hermes1D Release version: ${H1D_RELEASE}")
	message("\tBuild Hermes1D examples: ${H1D_WITH_EXAMPLES}")
	message("\tBuild Hermes1D benchmarks: ${H1D_WITH_BENCHMARKS}")
	message("\tBuild Hermes1D tests: ${H1D_WITH_TESTS}")
endif(WITH_H1D)

message("Build Hermes2D: ${WITH_H2D}")
if(WITH_H2D)
	message("\tBuild Hermes2D real version: ${H2D_REAL}")
	message("\tBuild Hermes2D complex version: ${H2D_COMPLEX}")
	message("\tBuild Hermes2D Debug version: ${H2D_DEBUG}")
	message("\tBuild Hermes2D Release version: ${H2D_RELEASE}")
	message("\tBuild Hermes2D with examples: ${H2D_WITH_EXAMPLES}")
	message("\tBuild Hermes2D with benchmarks: ${H2D_WITH_BENCHMARKS}")
	message("\tBuild Hermes2D with tests: ${H2D_WITH_TESTS}")
	message("\tBuild Hermes2D with tutorial: ${H2D_WITH_TUTORIAL}")
	message("\tBuild Hermes2D with GLUT: ${H2D_WITH_GLUT}")
	message("\tBuild Hermes2D with VIEWER_GUI: ${H2D_WITH_VIEWER_GUI}")
	message("\tBuild Hermes2D with util: ${H2D_WITH_UTIL}")
endif(WITH_H2D)

message("Build Hermes3D: ${WITH_H3D}")
if(WITH_H3D)
	message("\tBuild Hermes3D real version: ${H3D_REAL}")
	message("\tBuild Hermes3D complex version: ${H3D_COMPLEX}")
        message("\tBuild Hermes3D Debug version: ${H3D_DEBUG}")
   	message("\tBuild Hermes3D Release version: ${H3D_RELEASE}")
	message("\tBuild Hermes3D with Debug order: ${H3D_DEBUG_ORDER}")
	message("\tBuild Hermes3D with profiling: ${H3D_PROFILING}")
	message("\tBuild Hermes3D with Tetra: ${WITH_TETRA}")
	message("\tBuild Hermes3D with Hex: ${WITH_HEX}")
	message("\tBuild Hermes3D with Prism: ${WITH_PRISM}")
	message("\tBuild Hermes3D with examples: ${H3D_WITH_EXAMPLES}")
	message("\tBuild Hermes3D with benchmarks: ${H3D_WITH_BENCHMARKS}")
	message("\tBuild Hermes3D with tests: ${H3D_WITH_TESTS}")
	message("\tBuild Hermes3D with tools: ${H3D_WITH_TOOLS}")
	message("\tBuild Hermes3D with development tests: ${H3D_WITH_DEV_TESTS}")
endif(WITH_H3D)
message("\n")

#------------------------------------------------------------------------------
# Save the configuration.
#
# The following lines save the value of the configuration variables to a file
# that may be included by an external project dependent on Hermes. At this 
# moment, only variables needed for configuring our first modules are saved.
#
# NOTE: This command must be issued after all subdirectories with H1D,2D,3D
# have been processed, in order to include their definitions (alternatively, one
# could append these defs. to the config file from within the respective 
# CMakeLists.
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FILE(WRITE config_real.cmake 
"# HERMES-REAL CMAKE CONFIGURATION FILE\n
SET(WITH_UMFPACK        ${WITH_UMFPACK})
SET(WITH_SUPERLU        ${WITH_SUPERLU})
SET(SUPERLU_MT          ${SUPERLU_MT})
SET(WITH_PETSC          ${WITH_PETSC})
SET(WITH_MUMPS          ${WITH_MUMPS})
SET(WITH_TRILINOS       ${WITH_TRILINOS})
SET(WITH_EXODUSII       ${WITH_EXODUSII})
SET(WITH_HDF5           ${WITH_HDF5})
SET(WITH_OPENMP         ${WITH_OPENMP})

SET(HDF5_LIBRARY        ${HDF5_LIBRARY})
SET(UMFPACK_LIBRARIES   ${UMFPACK_LIBRARIES})
SET(TRILINOS_LIBRARIES  ${TRILINOS_LIBRARIES})
SET(PETSC_LIBRARIES     ${PETSC_REAL_LIBRARIES})
SET(MUMPS_LIBRARIES     ${MUMPS_REAL_LIBRARIES})
SET(SUPERLU_LIBRARY     ${SUPERLU_LIBRARY})
SET(LAPACK_LIBRARIES    ${LAPACK_LIBRARIES})
SET(PTHREAD_LIBRARY     ${PTHREAD_LIBRARY})
SET(CLAPACK_LIBRARY     ${CLAPACK_LIBRARY})
SET(BLAS_LIBRARY        ${BLAS_LIBRARY})
SET(F2C_LIBRARY         ${F2C_LIBRARY})
SET(ADDITIONAL_LIBS     ${ADDITIONAL_LIBS})
SET(PYTHON_LIBRARY      ${PYTHON_LIBRARY})
SET(EXODUSII_LIBRARIES  ${EXODUSII_LIBRARIES})
SET(WITH_GLUT           ${H2D_WITH_GLUT})
SET(GLUT_LIBRARY        ${GLUT_LIBRARY})
SET(GLEW_LIBRARY        ${GLEW_LIBRARY})

SET(UMFPACK_INCLUDE_DIR  ${UMFPACK_INCLUDE_DIR})
LIST(APPEND PETSC_INCLUDE_DIR ${COMMON_PETSC_INCLUDE_DIRS} ${PETSC_REAL_INCLUDE_DIRS})
SET(MUMPS_INCLUDE_DIR    ${MUMPS_INCLUDE_DIR})
SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR})
SET(SUPERLU_INCLUDE_DIR  ${SUPERLU_INCLUDE_DIR})
SET(HDF5_INCLUDE_DIR     ${HDF5_INCLUDE_DIR})
SET(EXODUSII_INCLUDE_DIR ${EXODUSII_INCLUDE_DIR})
SET(CLAPACK_INCLUDE_DIRS ${CLAPACK_INCLUDE_DIRS})
SET(NUMPY_INCLUDE_PATH   ${NUMPY_INCLUDE_PATH})
SET(PYTHON_INCLUDE_PATH  ${PYTHON_INCLUDE_PATH})
SET(DEFAULT_PYTHON_INSTALL_PATH  ${DEFAULT_PYTHON_INSTALL_PATH})

SET(HAVE_VASPRINTF      ${HAVE_VASPRINTF})
SET(HAVE_STRCASECMP     ${HAVE_STRCASECMP})
SET(HAVE_TEUCHOS_LINK   ${HAVE_TEUCHOS_LINK})
SET(HAVE_TEUCHOS_BFD    ${HAVE_TEUCHOS_BFD})
SET(HAVE_EXECINFO       ${HAVE_EXECINFO})
SET(HAVE_FMEMOPEN       ${HAVE_FMEMOPEN})
SET(HAVE_LOG2           ${HAVE_LOG2})
")

install( FILES config_real.cmake DESTINATION include/hermes )   
