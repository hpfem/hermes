set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/hermes_common/cmake/cxx_flag_overrides.cmake)
project(hermes)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
#
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/hermes_common)

# Testing.
set(WITH_TESTS                YES)

# Optional parts of the library:
set(WITH_H1D                  YES)
	# Release and debug versions.
	set(H1D_DEBUG             YES)
	set(H1D_RELEASE           NO)
	# Optional parts of the library.
	set(H1D_WITH_EXAMPLES     YES)
	set(H1D_WITH_BENCHMARKS   YES)
	set(H1D_WITH_TESTS        YES)

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
		
set(WITH_H3D                  YES)
	# Real/complex version of the library.
	set(H3D_REAL              YES)
	set(H3D_COMPLEX           YES)
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
	

# Logo will be shown:
set(REPORT_WITH_LOGO        YES)

#### Solvers ###
# Standard UMFPACK.
set(WITH_UMFPACK			YES)

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
set(PETSC_ROOT              /opt/petsc/petsc-3.1-p4)
set(PETSC_ARCH              petsc-arch-without-scalar-type)

# If PARDISO is enabled, WITH_OPENMP must be set as well (and perhaps some other libraries
# added to ADDITIONAL_LIBS according to PARDISO user guide).
#
set(WITH_PARDISO            NO)
set(PARDISO_LIB             /opt/pardiso/libpardiso400_GNU430_AMD_IA64.so)

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
set(WITH_TRILINOS			NO)


### Mesh format ###
option(WITH_EXODUSII "Enable support for EXODUSII mesh format" NO)
option(WITH_HDF5     "Enable support for HDF5" NO)


### Others ###
# Parallel execution (tells the linker to use parallel versions of the selected solvers, if available):
#
# Enable OpenMP for shared memory processing. Note that this is required when PARDISO is enabled.
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
set(WITH_PYTHON             NO)

# Allow to override the default values in CMake.vars:
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
  include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# Find necessary packages.
find_package(PTHREAD REQUIRED)

# Solvers.
if(MSVC)
	if(WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
		find_package(CLAPACK REQUIRED)
		include_directories(${CLAPACK_INCLUDE_DIRS})
	endif(WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
else(MSVC)
  # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
  # FindLAPACK is Fortran-only so Fortran must be enabled as well.
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    enable_language(Fortran)
    find_package(LAPACK REQUIRED)
    add_definitions(-DWITH_BLAS)
  endif(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
endif(MSVC)  

if(WITH_PETSC)
	find_package(PETSC REQUIRED)
	# For some header files, PETSc uses the same name for both real and complex versions (they only reside in 
	# different directories), so we can't globally set the architecture-specific include directory here, but 
	# rather on a per-target basis.
	include_directories(${COMMON_PETSC_INCLUDE_DIRS})
endif(WITH_PETSC)

if(WITH_UMFPACK)
	find_package(UMFPACK REQUIRED)
	include_directories(${UMFPACK_INCLUDE_DIRS})
endif(WITH_UMFPACK)

if(WITH_PARDISO)
	find_package(PARDISO REQUIRED)
endif(WITH_PARDISO)

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

add_subdirectory(${HERMES_COMMON_DIR})

include(hermes_common/CMakeCommon.txt)
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

if(WITH_TESTS)
  enable_testing()
endif(WITH_TESTS)

message("\nConfiguration results")
message("---------------------")
message("Build with Python: ${WITH_PYTHON}")
message("Build with UMFPACK: ${WITH_UMFPACK}")
message("Build with PETSC: ${WITH_PETSC}")
message("Build with PARDISO: ${WITH_PARDISO}")
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
message("Build Hermes1D: ${WITH_H1D}")
if(WITH_H1D)
	message("\tBuild Hermes1D Debug version: ${H1D_DEBUG}")
	message("\tBuild Hermes1D Release version: ${H1D_RELEASE}")
		# Optional parts of the library.
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
