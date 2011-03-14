cmake_minimum_required(VERSION 2.6)

project(hermes3d)

IF(H3D_DEBUG)
	set(DEBUG YES)
ELSE(H3D_DEBUG)
	set(DEBUG NO)
ENDIF(H3D_DEBUG)

IF(H3D_RELEASE)
  set(RELEASE YES)
  message(WARNING "As of now, Hermes3D is just an experimental package. Internal
                   compiler errors have been observed on GCC-4.4.1 -- if this is
                   your case, please set H3D_RELEASE to NO.") 
ELSE(H3D_RELEASE)
  set(RELEASE NO)
ENDIF(H3D_RELEASE)

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

#
# Version configuration.
#

if (NOT DEBUG AND NOT RELEASE)
  message(FATAL_ERROR "WITH_H3D is set to YES, but no versions are specified. 
                       Re-run with DEBUG and/or RELEASE set to YES.")
endif (NOT DEBUG AND NOT RELEASE)

if(H3D_REAL)
    list(APPEND VERS "real")
endif(H3D_REAL)
if(H3D_COMPLEX)
    list(APPEND VERS "complex")
endif(H3D_COMPLEX)

list(LENGTH VERS VERS_LEN)
if(NOT(VERS_LEN GREATER 0))
  message(FATAL_ERROR "WITH_H3D is set to YES, but no versions are specified. 
                       Re-run with REAL and/or COMPLEX set to YES.")
endif(NOT(VERS_LEN GREATER 0))

if(MSVC)
    if(H3D_REAL)
        set(HERMES3D_REAL "hermes3d-real")              # Base name of the library.
        set(HERMES3D_REAL_RELEASE ${HERMES3D_REAL})     # Name of the release version.
        set(HERMES3D_REAL_DEBUG ${HERMES3D_REAL}-debug) # Name of the debug version.
    endif(H3D_REAL)
    if(H3D_COMPLEX)
        set(HERMES3D_CPLX "hermes3d-cplx")
        set(HERMES3D_CPLX_RELEASE ${HERMES3D_CPLX})
        set(HERMES3D_CPLX_DEBUG ${HERMES3D_CPLX}-debug)
    endif(H3D_COMPLEX)
else(MSVC)
  if(DEBUG)
    if(H3D_REAL)
      set(HERMES3D_REAL_DEBUG "hermes3d-real-debug")
      set(HERMES3D_REAL ${HERMES3D_REAL_DEBUG})
    endif(H3D_REAL)
    if(H3D_COMPLEX)
      set(HERMES3D_CPLX_DEBUG "hermes3d-cplx-debug")
      set(HERMES3D_CPLX ${HERMES3D_CPLX_DEBUG})
    endif(H3D_COMPLEX)
    
    # set(CMAKE_BUILD_TYPE Debug) # This does not work with multi-version configurations,
                                  # because it is linked to the project (i.e. we would have
                                  # to have separate projects for Debug and Release)
  endif(DEBUG)

  if(RELEASE)
    if(H3D_REAL)
      set(HERMES3D_REAL_RELEASE "hermes3d-real")
      set(HERMES3D_REAL ${HERMES3D_REAL_RELEASE})
    endif(H3D_REAL)
    if(H3D_COMPLEX)
      set(HERMES3D_CPLX_RELEASE "hermes3d-cplx")
      set(HERMES3D_CPLX ${HERMES3D_CPLX_RELEASE})
    endif(H3D_COMPLEX)

    #set(CMAKE_BUILD_TYPE Release)
  endif(RELEASE)
endif(MSVC)

#
# Optional libraries and settings for H3D.
#

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

# H3D needs BLAS and LAPACK/CLAPACK. 
add_definitions(-DWITH_BLAS)

if(MSVC)
	if(NOT CLAPACK_FOUND)
		find_package(CLAPACK REQUIRED)
		include_directories(${CLAPACK_INCLUDE_DIRS})
	endif(NOT CLAPACK_FOUND)
	
	# If the library was found, prevent searching for it again in subsequent 
	# invocations of cmake .
	set(CLAPACK_FOUND   YES CACHE STRING "Lapack/BLAS library.")
else(MSVC)
  if(NOT LAPACK_FOUND)
    enable_language(Fortran)
    find_package(LAPACK REQUIRED)
    
  	# If the library was found, prevent searching for it again in subsequent 
	# invocations of cmake .
	set(LAPACK_FOUND    YES CACHE STRING "Lapack/BLAS library.")
  endif(NOT LAPACK_FOUND)
endif(MSVC)  

# Directory where the output will be written.
set(OUTPUT_DIR "" CACHE STRING "Directory were the output will be written")

#
# Python.
#

# We require Python for the mesh parser (it also needs NumPy)
if(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)
  find_package(Python REQUIRED)
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)    

if(WITH_PYTHON)
  #find_package(PythonLibs REQUIRED)
  set(CYTHON_INCLUDE_DIRECTORIES -I.)
  set(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${CMAKE_SOURCE_DIR}/hermes_common)
  add_subdirectory(python)
endif(WITH_PYTHON)

#
# Subdirectories.
#

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

#
# Documentation.
#

# TODO: Find out why this section has been commented out.

# Doxygen related
#set(DOXYGEN_BINARY doxygen)
#set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

#set(SPHINX_BUILD_BIN sphinx-build)

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


#
# Packaging.
#

# TODO: Test the packaging system and update it if neccessary.

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
