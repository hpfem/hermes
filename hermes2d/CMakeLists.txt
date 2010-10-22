set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/../hermes_common/cmake/cxx_flag_overrides.cmake)
project(hermes2d)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
#
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/../hermes_common)

# Build with debugging symbols:
set(DEBUG                   YES)
set(RELEASE                 NO)

# Real/complex version of the library (you can build either one or both):
set(H2D_REAL                YES)
set(H2D_COMPLEX             YES)

# Optional parts of the library:
set(WITH_EXAMPLES           YES)
set(WITH_BENCHMARKS         YES)
set(WITH_TUTORIAL           YES)
set(WITH_TESTS              YES)

# Optional functionality:
set(WITH_PYTHON             NO)
set(WITH_GLUT               YES)
set(WITH_VIEWER_GUI         NO)
set(WITH_UTIL               YES)

#
# Select which solvers will be available:
#
set(WITH_UMFPACK            YES)
set(WITH_TRILINOS           NO)

# If PETSc is enabled, PETSC_ROOT and PETSC_ARCH must be set as well;
# if PETSc is disabled, PETSC_ROOT and PETSC_ARCH will be ignored.
#
set(WITH_PETSC              NO)
set(PETSC_ROOT              /opt/petsc/petsc-3.1-p4)
set(PETSC_ARCH              petsc-arch-without-scalar-type)

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

# If PARDISO is enabled, WITH_OPENMP must be set as well (and perhaps some other libraries
# added to ADDITIONAL_LIBS according to PARDISO user guide).
#
set(WITH_PARDISO            NO)
set(PARDISO_LIB             /opt/pardiso/libpardiso400_GNU430_AMD_IA64.so)

#
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


# Additional data formats:
set(WITH_EXODUSII           NO)
set(WITH_HDF5               NO)

# Additional libraries required by some of the above:
# set(ADDITIONAL_LIBS       -lgfortran -lm)

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

# Allow to override the default values in CMake.vars:
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
  include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# Set overrides:
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

find_package(PTHREAD REQUIRED)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

if(WITH_GLUT)
  find_package(GLUT REQUIRED)
  find_package(GLEW REQUIRED)
else(WITH_GLUT)
  add_definitions(-DNOGLUT)
endif(WITH_GLUT)

if(WITH_VIEWER_GUI)
  find_package(ANTTWEAKBAR REQUIRED)
  include_directories(${ANTTWEAKBAR_INCLUDE_DIR})
endif(WITH_VIEWER_GUI)

if(WITH_EXODUSII)
  find_package(EXODUSII REQUIRED)
  include_directories(${EXODUSII_INCLUDE_DIR})
endif(WITH_EXODUSII)

if(WITH_HDF5)
	find_package(HDF5 REQUIRED)
	include_directories(${HDF5_INCLUDE_DIR})
endif(WITH_HDF5)

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

if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif(WITH_OPENMP)

if(NOT MSVC)
  # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
  # FindLAPACK is Fortran-only so Fortran must be enabled as well.
  if(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS)
    enable_language(Fortran)
    find_package(LAPACK REQUIRED)
    add_definitions(-DWITH_BLAS)
  endif(WITH_UMFPACK OR WITH_PETSC OR WITH_PARDISO OR WITH_MUMPS OR WITH_TRILINOS)
endif(NOT MSVC)  
 
# Solvers:
if(WITH_UMFPACK)
  find_package(UMFPACK REQUIRED)
  include_directories(${UMFPACK_INCLUDE_DIRS})
endif(WITH_UMFPACK)

if(WITH_TRILINOS)
  find_package(TRILINOS REQUIRED)
  include_directories(${TRILINOS_INCLUDE_DIR})
endif(WITH_TRILINOS)

if(WITH_MUMPS)
  find_package(MUMPS REQUIRED)
  if(WITH_MPI)
    find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS 
  endif(WITH_MPI)
  include_directories(${MUMPS_INCLUDE_DIR})
endif(WITH_MUMPS)

if(WITH_PETSC)
  find_package(PETSC REQUIRED)
  # For some header files, PETSc uses the same name for both real and complex versions (they only reside in 
  # different directories), so we can't globally set the architecture-specific include directory here, but 
  # rather on a per-target basis.
  include_directories(${COMMON_PETSC_INCLUDE_DIRS})
endif(WITH_PETSC)

if(WITH_PARDISO)
  find_package(PARDISO REQUIRED) 
endif(WITH_PARDISO)

include_directories(${HERMES_COMMON_DIR})

include(CheckIncludeFile)
CHECK_INCLUDE_FILE("link.h" HAVE_TEUCHOS_LINK)
CHECK_INCLUDE_FILE("bfd.h" HAVE_TEUCHOS_BFD)
CHECK_INCLUDE_FILE("execinfo.h" HAVE_EXECINFO)
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(HAVE_CXXABI yes)
ENDIF()
include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(vasprintf HAVE_VASPRINTF)

IF (HAVE_EXECINFO AND HAVE_CXXABI AND HAVE_VASPRINTF AND HAVE_TEUCHOS_LINK AND HAVE_TEUCHOS_BFD AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(HAVE_TEUCHOS_STACKTRACE yes)
ENDIF()

#add_subdirectory(hermes_common)
add_subdirectory(src)

if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(WITH_TESTS)

if(WITH_PYTHON)
  add_subdirectory(python)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
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

add_custom_target(test COMMAND ctest)
add_custom_target(test-quick COMMAND ctest -LE slow -j9)

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
message("Build with UMFPACK: ${WITH_UMFPACK}")
message("Build with PETSC: ${WITH_PETSC}")
message("Build with PARDISO: ${WITH_PARDISO}")
message("Build with MUMPS: ${WITH_MUMPS}")
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
