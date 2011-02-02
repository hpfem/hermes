project(hermes3d)

cmake_minimum_required(VERSION 2.6)
# H3D only supports Debug version.
option(DEBUG         "Build debug version" YES)

# Directory where the output will be written.
set(OUTPUT_DIR "" CACHE STRING "Directory were the output will be written")

IF(H3D_WITH_EXAMPLES)
	set(WITH_EXAMPLES YES)
ELSE(H3D_WITH_EXAMPLES)
	set(WITH_EXAMPLES NO)
ENDIF(H3D_WITH_EXAMPLES)

IF(H3D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS YES)
ELSE(H3D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS NO)
ENDIF(H3D_WITH_BENCHMARKS)

IF(H3D_WITH_TESTS)
	set(WITH_TESTS YES)
ELSE(H3D_WITH_TESTS)
	set(WITH_TESTS NO)
ENDIF(H3D_WITH_TESTS)

IF(H3D_WITH_TOOLS)
	set(WITH_TOOLS YES)
ELSE(H3D_WITH_TOOLS)
	set(WITH_TOOLS NO)
ENDIF(H3D_WITH_TOOLS)

IF(H3D_WITH_DEV_TESTS)
	set(WITH_DEV_TESTS YES)
ELSE(H3D_WITH_DEV_TESTS)
	set(WITH_DEV_TESTS NO)
ENDIF(H3D_WITH_DEV_TESTS)

# Doxygen related
#set(DOXYGEN_BINARY doxygen)
#set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

#set(SPHINX_BUILD_BIN sphinx-build)

# need at least real or complex version
if(H3D_REAL)
	list(APPEND VERS "real")
endif(H3D_REAL)
if(H3D_COMPLEX)
	list(APPEND VERS "complex")
endif(H3D_COMPLEX)

list(LENGTH VERS VERS_LEN)
if(NOT(VERS_LEN GREATER 0))
	message(FATAL_ERROR "H3D: No versions specified. Re-run with REAL and/or COMPLEX set to YES")
endif(NOT(VERS_LEN GREATER 0))

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
if(NOT(ELEMTYPES_LEN GREATER 0))
	message(FATAL_ERROR "H3D: No element types selected. Re-run with WITH_TETRA, WITH_HEX or WITH_PRISM set to YES")
endif(NOT(ELEMTYPES_LEN GREATER 0))

# internals regarding build type
if(DEBUG)
	if(H3D_REAL)
		set(HERMES3D_REAL hermes3d-debug)
	endif(H3D_REAL)
	if(H3D_COMPLEX) 
		set(HERMES3D_CPLX hermes3d-cplx-debug)
	endif(H3D_COMPLEX)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	if(H3D_REAL)
		set(HERMES3D_REAL hermes3d)
	endif(H3D_REAL)
	if(H3D_COMPLEX)
		set(HERMES3D_CPLX hermes3d-cplx)
	endif(H3D_COMPLEX)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)

# H3D needs BLAS. However, if it was needed before, there is no need to look for it again.
add_definitions(-DWITH_BLAS)
if(MSVC)
	if(NOT(WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU))
		find_package(CLAPACK REQUIRED)
		include_directories(${CLAPACK_INCLUDE_DIRS})
	endif(NOT(WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU))
else(MSVC)
  if(NOT(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU))
    enable_language(Fortran)
    find_package(LAPACK REQUIRED)
  endif(NOT(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU))
endif(MSVC)  

# We require Python for the mesh parser (it also needs NumPy)
find_package(Python REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
set(PYTHON_INSTALL_PATH lib/python/site-packages)

# Python inclusion is different in H1D, H2D, H3D.
if(WITH_PYTHON)
	find_package(PythonLibs REQUIRED)
	find_package(Cython REQUIRED)
	add_subdirectory(python)
	include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(WITH_PYTHON)

# Subdirectories.
add_subdirectory(src)

if(WITH_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)
if(WITH_TOOLS)
	add_subdirectory(tools)
endif(WITH_TOOLS)

# examples
add_subdirectory(examples)
add_custom_target(examples)
add_custom_command(
	TARGET    examples
	COMMAND   make -C ${PROJECT_BINARY_DIR}/examples
)

# benchmarks
add_subdirectory(benchmarks)
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