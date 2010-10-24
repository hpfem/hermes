project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)

set(WITH_PYTHON yes)
set(DEBUG yes)
set(RELEASE no)
set(WITH_EXAMPLES yes)
set(WITH_TESTS yes)

option(WITH_UMFPACK  "Enable support for UMFPACK" YES)
option(WITH_PETSC    "Enable support for PETSc" NO)
option(WITH_PARDISO  "Enable support for PARDISO" NO)
option(WITH_MUMPS    "Enable support for MUMPS solver" NO)
option(WITH_TRILINOS "Enable support for Trilinos solvers" NO)

option(WITH_EXODUSII "Enable support for EXODUSII mesh format" NO)
option(WITH_HDF5     "Enable support for HDF5" NO)


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

# For use in the common code.
set(H1D_REAL YES)
#add_definitions(-DH1D_REAL)


if(DEBUG)
        set(HERMES_BIN hermes1d-debug)
endif(DEBUG)

if(RELEASE)
        set(HERMES_BIN hermes1d)
endif(RELEASE)

add_subdirectory(src)


find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

if(WITH_PYTHON)
    include(UseCython)

    # if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
    # try to guess using cmake/FindPythonLibs.cmake
    if (NOT PYTHON_INCLUDE_PATH)
        find_package(PythonLibs REQUIRED)
    endif (NOT PYTHON_INCLUDE_PATH)
    find_package(NumPy REQUIRED)

    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
    include_directories(src)
    add_subdirectory(hermes1d)
endif(WITH_PYTHON)
if(WITH_EXAMPLES)
    add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

# find necessary packages
find_package(JUDY REQUIRED)
find_package(PTHREAD REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
include_directories(${HERMES_COMMON_DIR})


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
