set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/../hermes_common/cmake/cxx_flag_overrides.cmake) # this is in H2D
project(hermes3d)

cmake_minimum_required(VERSION 2.6)

# default values

option(DEBUG         "Build debug version" YES)
set(DEBUG_ORDER NO)
set(PROFILING   NO)

# real/complex version of the library
option(H3D_REAL          "Build real (double precision) version of Hermes3D" YES)
option(H3D_COMPLEX       "Build complex version of Hermes3D" YES)
# supported element types
option(WITH_TETRA    "Enable tetrahedral elements" YES)
option(WITH_HEX      "Enable hexahedral elements" YES)
option(WITH_PRISM    "Enable prism elements" NO)

option(WITH_UMFPACK  "Enable support for UMFPACK" YES)
option(WITH_PETSC    "Enable support for PETSc" NO)
option(WITH_PARDISO  "Enable support for PARDISO" NO)
option(WITH_MUMPS    "Enable support for MUMPS solver" NO)
option(WITH_TRILINOS "Enable support for Trilinos solvers" NO)

option(WITH_EXODUSII "Enable support for EXODUSII mesh format" NO)
option(WITH_HDF5     "Enable support for HDF5" NO)

option(WITH_OPENMP   "Build with OpenMP support" NO)

option(WITH_EXAMPLES "Build examples" YES)
option(WITH_BENCHMARKS "Build benchmarks" YES)
option(WITH_TESTS    "Build tests" YES)
set(WITH_TOOLS  NO)

option(WITH_PYTHON   "Enable support for Python" NO)

set(OUTPUT_DIR "" CACHE STRING "Directory were the output will be written")

set(DEV_TESTS NO)

# Doxygen related
#set(DOXYGEN_BINARY doxygen)
#set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

#set(SPHINX_BUILD_BIN sphinx-build)

# allow to override default values via CMake.vars 
include(${PROJECT_BINARY_DIR}/CMake.vars OPTIONAL)

# sanity checks

# need at least real or complex version
if(H3D_REAL)
	list(APPEND VERS "real")
endif(H3D_REAL)
if(H3D_COMPLEX)
	list(APPEND VERS "complex")
endif(H3D_COMPLEX)

list(LENGTH VERS VERS_LEN)
if(VERS_LEN GREATER 0)
	message(STATUS "Will build versions: ${VERS}")
else(VERS_LEN GREATER 0)
	message(FATAL_ERROR "No versions specified. Re-run with REAL and/or COMPLEX set to YES")
endif(VERS_LEN GREATER 0)

if(WITH_TETRA)
	list(APPEND ELEMTYPES "tetra")
endif(WITH_TETRA)
if(WITH_HEX)
	list(APPEND ELEMTYPES "hex")
endif(WITH_HEX)
if(WITH_PRISM)
	list(APPEND ELEMTYPES "prism")
endif(WITH_PRISM)

list(LENGTH ELEMTYPES ELEMTYPES_LEN)
if(ELEMTYPES_LEN GREATER 0)
	message(STATUS "Will build with element types: ${ELEMTYPES}")
else(ELEMTYPES_LEN GREATER 0)
	message(FATAL_ERROR "No element types selected. Re-run with WITH_TETRA, WITH_HEX or WITH_PRISM set to YES")
endif(ELEMTYPES_LEN GREATER 0)

# set internal library-wide variables
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/../hermes_common)
set(MESHES_DIR ${PROJECT_SOURCE_DIR}/meshes)

# internals regarding build type
if(DEBUG)
	if(H3D_REAL)
		set(HERMES_REAL_BIN hermes3d-debug)
	endif(H3D_REAL)
	if(H3D_COMPLEX) 
		set(HERMES_CPLX_BIN hermes3d-cplx-debug)
	endif(H3D_COMPLEX)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	if(H3D_REAL)
		set(HERMES_REAL_BIN hermes3d)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		set(HERMES_CPLX_BIN hermes3d-cplx)
	endif(H3D_COMPLEX)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)

# find necessary packages
find_package(JUDY REQUIRED)
find_package(PTHREAD REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
add_definitions(-DWITH_BLAS)

if(NOT MSVC)
	enable_language(Fortran)	# need a fortran compiler for checking BLAS
	find_package(LAPACK REQUIRED) # also searches for BLAS
endif(NOT MSVC)

# solvers
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

# SMP
if(WITH_OPENMP)
	find_package(OpenMP REQUIRED)
	set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif(WITH_OPENMP)


# domain decomposition

if(WITH_METIS)
	find_package(METIS REQUIRED)
endif(WITH_METIS)


if(WITH_HDF5)
	find_package(HDF5 REQUIRED)
	include_directories(${HDF5_INCLUDE_DIR})
endif(WITH_HDF5)

if(WITH_EXODUSII)
	find_package(EXODUSII REQUIRED)
	include_directories(${EXODUSII_INCLUDE_DIR})
endif(WITH_EXODUSII)

# If using any package that requires MPI (e.g. parallel versions of MUMPS, PETSC)
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

if(WITH_GLUT)
  find_package(GLUT REQUIRED)
  find_package(GLEW REQUIRED)
else(WITH_GLUT)
  add_definitions(-DNOGLUT)
endif(WITH_GLUT)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${HERMES_COMMON_DIR})

include(CheckIncludeFile)
CHECK_INCLUDE_FILE("link.h" HAVE_TEUCHOS_LINK)
CHECK_INCLUDE_FILE("bfd.h" HAVE_TEUCHOS_BFD)

IF (HAVE_TEUCHOS_LINK AND HAVE_TEUCHOS_BFD AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(HAVE_TEUCHOS_STACKTRACE yes)
ENDIF()

#
add_subdirectory(src)
if(WITH_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)
if(WITH_TOOLS)
	add_subdirectory(tools)
endif(WITH_TOOLS)

if(WITH_PYTHON)
	find_package(PythonLibs REQUIRED)
	find_package(Cython REQUIRED)
	add_subdirectory(python)
	message(STATUS "Compiling with python support")
	include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(WITH_PYTHON)


# examples
add_subdirectory(examples EXCLUDE_FROM_ALL)
add_custom_target(examples)
add_custom_command(
	TARGET    examples
	COMMAND   make -C ${PROJECT_BINARY_DIR}/examples
)

# benchmarks
add_subdirectory(benchmarks EXCLUDE_FROM_ALL)
add_custom_target(benchmarks)
add_custom_command(
	TARGET    benchmarks
	COMMAND   make -C ${PROJECT_BINARY_DIR}/benchmarks
)

# doc
#add_subdirectory(doc)

#add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/xml/index.xml)
#add_custom_command(
#	TARGET    doc
#	COMMAND   ${SPHINX_BUILD_BIN} -N -b html ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR}/doc/html
#)

#add_custom_target(doc-tex DEPENDS ${PROJECT_BINARY_DIR}/doc/xml/index.xml)
#add_custom_command(
#	TARGET    doc-tex
#	COMMAND   ${SPHINX_BUILD_BIN} -N -b latex ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR}/doc/tex
#)

#add_custom_command(
#	OUTPUT    ${PROJECT_BINARY_DIR}/doc/xml/index.xml
#	COMMAND   ${DOXYGEN_BINARY} ${DOXYGEN_CONFIG_FILE}
#)


add_custom_target(test-quick
	COMMAND   /usr/bin/ctest -LE long
)

# -- package --------

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_PACKAGE_VENDOR "hp-FEM Group, University of Nevada, Reno")
set(CPACK_SYSTEM_NAME "linux")

# source
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "hermes3d-${VERSION}")
set(CPACK_SOURCE_IGNORE_FILES 
	"${PROJECT_SOURCE_DIR}/.hg/"
	"${PROJECT_SOURCE_DIR}/.settings/"
	"${PROJECT_SOURCE_DIR}/.hgignore"
	"${PROJECT_SOURCE_DIR}/.project"
	"${PROJECT_SOURCE_DIR}/.cproject"
	"${PROJECT_SOURCE_DIR}/build/"
	"${PROJECT_SOURCE_DIR}/CMake.vars$"
	"${PROJECT_SOURCE_DIR}/utils/"
)

include(CPack)

if(HAVE_TEUCHOS_STACKTRACE)
    message("-- Print Teuchos stacktrace on segfault: YES")
else(HAVE_TEUCHOS_STACKTRACE)
    message("-- Print Teuchos stacktrace on segfault: NO")
endif(HAVE_TEUCHOS_STACKTRACE)
