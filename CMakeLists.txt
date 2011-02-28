
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/hermes_common/cmake/cxx_flag_overrides.cmake)
project(hermes)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/hermes_common)
set(HERMES_COMMON_LIB_REAL  hermes_common-real)
set(HERMES_COMMON_LIB_CPLX  hermes_common-cplx)
#------------------------------------------------------------------------------
# Main configuration options.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### BEGIN DEFAULT CONFIGURATION ###

# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

# H1D:
set(WITH_H1D                  YES)
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
	set(H2D_WITH_MODULES      YES)
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

# Logo will be shown:
set(REPORT_WITH_LOGO        YES)

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
  set(PYTHON_INSTALL_PATH lib/python/site-packages)
endif(WITH_PYTHON)

# Solvers (TODO: This could maybe go into hermes_common/solver/CMakeLists.txt).
if(MSVC)
	if(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
	    if (NOT CLAPACK_FOUND)
    		find_package(CLAPACK REQUIRED)
	    	include_directories(${CLAPACK_INCLUDE_DIRS})
	    endif (NOT CLAPACK_FOUND)
	 
	    # If the library was found, prevent searching for it again in subsequent 
	    # invocations of cmake .
	    set(CLAPACK_FOUND    YES CACHE STRING "Lapack/BLAS library.")   
	endif(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
else(MSVC)
  # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
  # FindLAPACK is Fortran-only so Fortran must be enabled as well.
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    if (NOT LAPACK_FOUND)
        enable_language(Fortran)
        find_package(LAPACK REQUIRED)
    endif (NOT LAPACK_FOUND)
    
    add_definitions(-DWITH_BLAS)
    
    # If the library was found, prevent searching for it again in subsequent 
	# invocations of cmake .
	set(LAPACK_FOUND    YES CACHE STRING "Lapack/BLAS library.")
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
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # Do we have to add this in 
                                                # set_target_properties?
endif(MSVC)

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
message("---------------------")
message("Hermes common library:")
message("\tBuild Real version: ${HERMES_COMMON_REAL}")
message("\tBuild Complex version: ${HERMES_COMMON_COMPLEX}")
message("\tBuild Debug version: ${HERMES_COMMON_DEBUG}")
message("\tBuild Release version: ${HERMES_COMMON_RELEASE}")
message("\tBuild with tests: ${HERMES_COMMON_WITH_TESTS}")

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
	message("\tBuild Hermes2D examples: ${H2D_WITH_EXAMPLES}")
	message("\tBuild Hermes2D benchmarks: ${H2D_WITH_BENCHMARKS}")
	message("\tBuild Hermes2D tests: ${H2D_WITH_TESTS}")
	message("\tBuild Hermes2D modules: ${H2D_WITH_MODULES}")
	message("\tBuild Hermes2D tutorial: ${H2D_WITH_TUTORIAL}")
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
	message("\tBuild Hermes3D examples: ${H3D_WITH_EXAMPLES}")
	message("\tBuild Hermes3D with benchmarks: ${H3D_WITH_BENCHMARKS}")
	message("\tBuild Hermes3D with tests: ${H3D_WITH_TESTS}")
	message("\tBuild Hermes3D with tools: ${H3D_WITH_TOOLS}")
	message("\tBuild Hermes3D with development tests: ${H3D_WITH_DEV_TESTS}")
endif(WITH_H3D)
message("\n")
