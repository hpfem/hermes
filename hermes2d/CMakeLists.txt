set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cxx_flag_overrides.cmake)
project(hermes2d)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# build with debugging symbols
set(DEBUG YES)
set(RELEASE NO)

# real/complex version of the library (you can build either of both)
set(H2D_REAL     YES)
set(H2D_COMPLEX  YES)

# optional functionality
set(WITH_PYTHON     NO)
set(WITH_EXAMPLES   YES)
set(WITH_BENCHMARKS YES)
set(WITH_TUTORIAL   YES)
set(WITH_TESTS      YES)
set(WITH_GLUT       YES)
set(WITH_VIEWER_GUI NO)
set(WITH_UTIL       YES)
set(WITH_TRILINOS   NO)
set(WITH_EXODUSII   NO)

# reporting and logging
set(REPORT_WITH_LOGO YES) #logo will be shown
set(REPORT_TO_FILE  YES) #reports will be mirrored to a file
set(REPORT_RUNTIME_CONTROL NO) #reporting will be controlled by preprocessor directives. The runtime control is inteded just for Python wrappers.
set(REPORT_ALL     NO) #do not report all events, just selected ones
set(REPORT_WARNING YES) #warnings will be reported
set(REPORT_INTR_WARNING NO) #warnings about itegration points will NOT be reported
set(REPORT_INFO     YES) #info will not be reported
set(REPORT_VERBOSE  NO) #info details will not be reported
set(REPORT_TRACE    NO) #code execution tracing will not be reported
set(REPORT_TIME     NO) #time will not be measured and time measurement will not be reported
#set(REPORT_DEBUG    NO) #debug events will depend on version which is compiled

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# set overrides
if(WITH_PYTHON)
	set(REPORT_RUNTIME_CONTROL YES) #reuired by python wrappers
endif(WITH_PYTHON)

# end of configuration section
#---------------------------------------------------------

if(WITH_TRILINOS)
  enable_language(Fortran)
  find_package(LAPACK REQUIRED)     # trilinos needs LAPACK
  find_package(TRILINOS REQUIRED)
endif(WITH_TRILINOS)

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

if(NOT WITH_GLUT)
    add_definitions(-DNOGLUT)
endif(NOT WITH_GLUT)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(WITH_EXODUSII)
    find_package(EXODUSII REQUIRED)
endif(WITH_EXODUSII)

add_subdirectory(hermes_common)
add_subdirectory(src)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

if(WITH_PYTHON)
	add_subdirectory(python)
endif(WITH_PYTHON)

if(WITH_EXAMPLES)
	add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
	add_subdirectory(benchmarks)
endif(WITH_BENCHMARKS)

if(WITH_TUTORIAL)
	add_subdirectory(tutorial)
endif(WITH_TUTORIAL)

if(WITH_UTIL)
	add_subdirectory(util)
endif(WITH_UTIL)

add_custom_target(test
    COMMAND ctest
    )
add_custom_target(test-quick
    COMMAND ctest -LE slow -j9
    )

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
message("Build with tutorial: ${WITH_TUTORIAL}")
message("Build with GLUT: ${WITH_GLUT}")
message("Build with viewer GUI: ${WITH_VIEWER_GUI}")
message("Build with util: ${WITH_UTIL}")
message("Build with tests: ${WITH_TESTS}")
message("Build with TRILINOS: ${WITH_TRILINOS}")
message("---------------------")
message("Hermes2D logo: ${REPORT_WITH_LOGO}")
message("Mirror reports to a log file: ${REPORT_TO_FILE}")
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
