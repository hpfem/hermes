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
        set(HERMES_REAL_BIN hermes2d-real)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
        set(HERMES_CPLX_BIN hermes2d-cplx)
    endif(H2D_COMPLEX)
else(MSVC)
  if(DEBUG)
    if(H2D_REAL)
      set(HERMES_REAL_BIN hermes2d-real-debug)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES_CPLX_BIN hermes2d-cplx-debug)
    endif(H2D_COMPLEX)
  endif(DEBUG)

  if(RELEASE)
    if(H2D_REAL)
      set(HERMES_REAL_BIN hermes2d-real)
    endif(H2D_REAL)
    if(H2D_COMPLEX)
      set(HERMES_CPLX_BIN hermes2d-cplx)
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


message("\nConfiguration results (these options can be set on the commandline)")
message("---------------------")
message("Versions to build: ${VERS}")
message("Build with Python: ${WITH_PYTHON}")
message("Build with examples: ${WITH_EXAMPLES}")
message("Build with benchmarks: ${WITH_BENCHMARKS}")
message("Build with modules: ${WITH_MODULES}")
message("Build with tutorial: ${WITH_TUTORIAL}")
message("Build with GLUT: ${H2D_WITH_GLUT}")
message("Build with viewer GUI: ${H2D_WITH_VIEWER_GUI}")
message("Build with util: ${H2D_WITH_UTIL}")
message("Build with tests: ${WITH_TESTS}")
message("Build with UMFPACK: ${WITH_UMFPACK}")
message("Build with PETSC: ${WITH_PETSC}")
message("Build with PARDISO: ${WITH_PARDISO}")
message("Build with MUMPS: ${WITH_MUMPS}")
message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
message("Build with TRILINOS: ${WITH_TRILINOS}")
message("Build with MPI: ${WITH_MPI}")
message("Build with OPENMP: ${WITH_OPENMP}")
message("---------------------")
message("Hermes2D logo: ${REPORT_WITH_LOGO}")
message("Mirror reports to a log file: ${REPORT_TO_FILE}")
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

message("\n")
