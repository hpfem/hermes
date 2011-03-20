cmake_minimum_required(VERSION 2.6)

project(hermes2d)

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)
#set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/../hermes_common)

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

IF(H2D_WITH_EXAMPLES)
	set(WITH_EXAMPLES YES)
ELSE(H2D_WITH_EXAMPLES)
	set(WITH_EXAMPLES NO)
ENDIF(H2D_WITH_EXAMPLES)

IF(H2D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS YES)
ELSE(H2D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS NO)
ENDIF(H2D_WITH_BENCHMARKS)

IF(H2D_WITH_TESTS)
	set(WITH_TESTS YES)
ELSE(H2D_WITH_TESTS)
	set(WITH_TESTS NO)
ENDIF(H2D_WITH_TESTS)

IF(H2D_WITH_TUTORIAL)
	set(WITH_TUTORIAL YES)
ELSE(H2D_WITH_TUTORIAL)
	set(WITH_TUTORIAL NO)
ENDIF(H2D_WITH_TUTORIAL)

#
# Version configuration.
#

if (NOT DEBUG AND NOT RELEASE)
  message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified. 
                       Re-run with DEBUG and/or RELEASE set to YES.")
endif (NOT DEBUG AND NOT RELEASE)

if(H2D_REAL)
    list(APPEND VERS "real")
endif(H2D_REAL)
if(H2D_COMPLEX)
    list(APPEND VERS "complex")
endif(H2D_COMPLEX)

list(LENGTH VERS VERS_LEN)
if(NOT(VERS_LEN GREATER 0))
  message(FATAL_ERROR "WITH_H2D is set to YES, but no versions are specified. 
                       Re-run with REAL and/or COMPLEX set to YES")
endif(NOT(VERS_LEN GREATER 0))

if(MSVC)
    if(H2D_REAL)
        set(HERMES2D_REAL "hermes2d-real")              # Base name of the library.
        set(HERMES2D_REAL_RELEASE ${HERMES2D_REAL})     # Name of the release version.
        set(HERMES2D_REAL_DEBUG ${HERMES2D_REAL}-debug) # Name of the debug version.
    endif(H2D_REAL)
    if(H2D_COMPLEX)
        set(HERMES2D_CPLX "hermes2d-cplx")
        set(HERMES2D_CPLX_RELEASE ${HERMES2D_CPLX})
        set(HERMES2D_CPLX_DEBUG ${HERMES2D_CPLX}-debug)
    endif(H2D_COMPLEX)
else(MSVC)
  if(DEBUG)
    if(H2D_REAL)
      set(HERMES2D_REAL_DEBUG "hermes2d-real-debug")
      set(HERMES2D_REAL ${HERMES2D_REAL_DEBUG})
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES2D_CPLX_DEBUG "hermes2d-cplx-debug")
      set(HERMES2D_CPLX ${HERMES2D_CPLX_DEBUG})
    endif(H2D_COMPLEX)
    
    # set(CMAKE_BUILD_TYPE Debug) # This does not work with multi-version configurations,
                                  # because it is linked to the project (i.e. we would have
                                  # to have separate projects for Debug and Release)
  endif(DEBUG)

  if(RELEASE)
    if(H2D_REAL)
      set(HERMES2D_REAL_RELEASE "hermes2d-real")
      set(HERMES2D_REAL ${HERMES2D_REAL_RELEASE})
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES2D_CPLX_RELEASE "hermes2d-cplx")
      set(HERMES2D_CPLX ${HERMES2D_CPLX_RELEASE})
    endif(H2D_COMPLEX)

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

#
# Python
#

# We require Python for the mesh parser (it also needs NumPy)
if(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)
  find_package(Python REQUIRED)
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)    

if(WITH_PYTHON)
  set(REPORT_RUNTIME_CONTROL YES) # Required by python wrappers.
  set(CYTHON_INCLUDE_DIRECTORIES -I.)
  set(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${CMAKE_SOURCE_DIR}/hermes_common)
  add_subdirectory(hermes2d)
endif(WITH_PYTHON)

#
# Subdirectories.
#

add_subdirectory(src)

if(WITH_TESTS)
  enable_testing()  
  add_subdirectory(tests)
endif(WITH_TESTS)

if(WITH_EXAMPLES)
  add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
  add_subdirectory(benchmarks-general)
  add_subdirectory(benchmarks-nist)
endif(WITH_BENCHMARKS)

if(WITH_TUTORIAL)
  add_subdirectory(tutorial)
endif(WITH_TUTORIAL)

if(H2D_WITH_UTIL)
  add_subdirectory(util)
endif(H2D_WITH_UTIL)

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
