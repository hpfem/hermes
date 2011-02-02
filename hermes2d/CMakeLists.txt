project(hermes2d)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
#
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/../hermes_common)

IF(H2D_DEBUG)
	set(DEBUG YES)
ELSE(H2D_DEBUG)
	set(DEBUG NO)
ENDIF(H2D_DEBUG)

IF(H2D_RELEASE)
	set(RELEASE NO)
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

IF(H2D_WITH_MODULES)
	set(WITH_MODULES YES)
ELSE(H2D_WITH_MODULES)
	set(WITH_MODULES NO)
ENDIF(H2D_WITH_MODULES)

IF(H2D_WITH_TUTORIAL)
	set(WITH_TUTORIAL YES)
ELSE(H2D_WITH_TUTORIAL)
	set(WITH_TUTORIAL NO)
ENDIF(H2D_WITH_TUTORIAL)

# Python inclusion is different in H1D, H2D, H3D.
if(WITH_PYTHON)
  set(REPORT_RUNTIME_CONTROL YES) #reuired by python wrappers
endif(WITH_PYTHON)

#
# end of configuration section
#---------------------------------------------------------
if(H2D_REAL)
    list(APPEND VERS "real")
endif(H2D_REAL)
if(H2D_COMPLEX)
    list(APPEND VERS "complex")
endif(H2D_COMPLEX)

if(MSVC)
    if(H2D_REAL)
        set(HERMES2D_REAL hermes2d-real)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
        set(HERMES2D_CPLX hermes2d-cplx)
    endif(H2D_COMPLEX)
else(MSVC)
  if(DEBUG)
    if(H2D_REAL)
      set(HERMES2D_REAL hermes2d-real-debug)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES2D_CPLX hermes2d-cplx-debug)
    endif(H2D_COMPLEX)
  endif(DEBUG)

  if(RELEASE)
    if(H2D_REAL)
      set(HERMES2D_REAL hermes2d-real)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES2D_CPLX hermes2d-cplx)
    endif(H2D_COMPLEX)
  endif(RELEASE)
endif(MSVC)

if(H2D_WITH_GLUT)
  find_package(GLUT REQUIRED)
  find_package(GLEW REQUIRED)
else(H2D_WITH_GLUT)
  add_definitions(-DNOGLUT)
endif(H2D_WITH_GLUT)

if(H2D_WITH_VIEWER_GUI)
  find_package(ANTTWEAKBAR REQUIRED)
  include_directories(${ANTTWEAKBAR_INCLUDE_DIR})
endif(H2D_WITH_VIEWER_GUI)

# We require Python for the mesh parser (it also needs NumPy)
find_package(Python REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
set(PYTHON_INSTALL_PATH lib/python/site-packages)

if(WITH_PYTHON)
  find_package(Cython REQUIRED)
  add_subdirectory(hermes2d)
endif(WITH_PYTHON)

# Subdirectories.
add_subdirectory(src)

if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(WITH_TESTS)


if(WITH_EXAMPLES)
  add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
  add_subdirectory(benchmarks)
endif(WITH_BENCHMARKS)

if(WITH_MODULES)
  add_subdirectory(modules)
endif(WITH_MODULES)

if(WITH_TUTORIAL)
  add_subdirectory(tutorial)
endif(WITH_TUTORIAL)

if(H2D_WITH_UTIL)
  add_subdirectory(util)
endif(H2D_WITH_UTIL)

# Documentation
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
