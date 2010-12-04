project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)

set(WITH_PYTHON no)
set(DEBUG yes)
set(RELEASE no)
set(WITH_EXAMPLES yes)
set(WITH_BENCHMARKS yes)
set(WITH_TESTS yes)

option(WITH_UMFPACK  "Enable support for UMFPACK" YES)
option(WITH_SUPERLU  "Enable support for SUPERLU" NO)
option(WITH_PETSC    "Enable support for PETSc" NO)
option(WITH_PARDISO  "Enable support for PARDISO" NO)
option(WITH_MUMPS    "Enable support for MUMPS solver" NO)
option(WITH_TRILINOS "Enable support for Trilinos solvers" NO)

option(WITH_EXODUSII "Enable support for EXODUSII mesh format" NO)
option(WITH_HDF5     "Enable support for HDF5" NO)

option(WITH_OPENMP   "Build with OpenMP support" NO)

# Reporting and logging:
set(REPORT_WITH_LOGO        YES)  #logo will be shown
set(REPORT_TO_FILE          YES)  #reports will be mirrored to a file
set(REPORT_RUNTIME_CONTROL  NO)   #reporting will be controlled by preprocessor directives. The runtime control is inteded just for Python wrappers.
set(REPORT_ALL              NO)   #do not report all events, just selected ones
set(REPORT_WARNING          YES)  #warnings will be reported
set(REPORT_INTR_WARNING     NO)   #warnings about itegration points will NOT be reported
set(REPORT_INFO             YES)  #info will be reported
set(REPORT_VERBOSE          NO)   #info details will not be reported
set(REPORT_TRACE            NO)   #code execution tracing will not be reported
set(REPORT_TIME             NO)   #time will not be measured and time measurement will not be reported
#set(REPORT_DEBUG           NO)   #debug events will depend on version which is compiled


# Doxygen related
set(DOXYGEN_BINARY doxygen)
set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

set(SPHINX_BUILD_BIN sphinx-build)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# Location of the common library.
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/../hermes_common)

include(../hermes_common/CMakeCommon.txt)

# For use in the common code.
set(H1D_REAL YES)
add_definitions(-DH1D_REAL)


if(DEBUG)
        set(HERMES_BIN hermes1d-debug)
endif(DEBUG)

if(RELEASE)
        set(HERMES_BIN hermes1d)
endif(RELEASE)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

if(WITH_PYTHON)
    # if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
    # try to guess using cmake/FindPythonLibs.cmake
    if (NOT PYTHON_INCLUDE_PATH)
        find_package(PythonLibs REQUIRED)
    endif (NOT PYTHON_INCLUDE_PATH)
    find_package(NumPy REQUIRED)

    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(WITH_PYTHON)

if(WITH_EXAMPLES)
    add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
    add_subdirectory(benchmarks)
endif(WITH_BENCHMARKS)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

# find necessary packages
find_package(JUDY REQUIRED)
find_package(PTHREAD REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
include_directories(${HERMES_COMMON_DIR})

if(NOT MSVC)
  # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
  # FindLAPACK is Fortran-only so Fortran must be enabled as well.
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    enable_language(Fortran)
    find_package(LAPACK REQUIRED)
    add_definitions(-DWITH_BLAS)
  endif(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
endif(NOT MSVC)

# Solvers:
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

# SMP
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif(WITH_OPENMP)

add_subdirectory(src)

add_subdirectory(doc.old)

add_custom_target(doc.old DEPENDS ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml)
add_custom_command(
	TARGET    doc.old
	COMMAND   ${SPHINX_BUILD_BIN} -N -b html ${PROJECT_SOURCE_DIR}/doc.old ${PROJECT_BINARY_DIR}/doc.old/html
)

add_custom_target(doc-tex DEPENDS ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml)
add_custom_command(
	TARGET    doc-tex
	COMMAND   ${SPHINX_BUILD_BIN} -N -b latex ${PROJECT_SOURCE_DIR}/doc.old ${PROJECT_BINARY_DIR}/doc.old/tex
)

add_custom_command(
	OUTPUT    ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml
	COMMAND   ${DOXYGEN_BINARY} ${DOXYGEN_CONFIG_FILE}
)


add_custom_target(test-quick
	COMMAND   /usr/bin/ctest -E \"hnnd|qorder|lobatto|hex-h1-1face\"
)

message("\nConfiguration results:")
message("----------------------")
message("Build with Python: ${WITH_PYTHON}")
message("Build with debug: ${DEBUG}")
message("Build with release: ${RELEASE}")
message("Build with tests: ${WITH_TESTS}")
message("\n")
