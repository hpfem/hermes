cmake_minimum_required(VERSION 2.6)

project(hermes_common)

IF(HERMES_COMMON_WITH_TESTS)
        set(WITH_TESTS YES)
ELSE(HERMES_COMMON_WITH_TESTS)
        set(WITH_TESTS NO)
ENDIF(HERMES_COMMON_WITH_TESTS)

# Check for stacktrace and math library.
include(CMakeCommon.txt)

# For selecting real or complex libraries and include directories for MUMPS and
# PETSC.
include(PickRealOrCplxLibs)

# Additional building steps required by MSVC.
include(BuildAndInstallScripts)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Teuchos RCP
if (NOT WITH_TRILINOS)
  set( TEUCHOS_SOURCES
    third_party_codes/trilinos-teuchos/Teuchos_TypeNameTraits.cpp
    third_party_codes/trilinos-teuchos/Teuchos_dyn_cast.cpp
    third_party_codes/trilinos-teuchos/Teuchos_Ptr.cpp
    third_party_codes/trilinos-teuchos/Teuchos_RCPNode.cpp
    third_party_codes/trilinos-teuchos/Teuchos_TestForException.cpp
    third_party_codes/trilinos-teuchos/Teuchos_TypeNameTraits.cpp
  )
endif (NOT WITH_TRILINOS)

#
# Source files for the Hermes_common library.
#
set(SRC
  hermes_logging.cpp
  common_time_period.cpp
  callstack.cpp
  error.cpp
  utils.cpp
  matrix.cpp
  tables.cpp
  qsort.cpp
  third_party_codes/trilinos-teuchos/Teuchos_stacktrace.cpp
  solver/nox.cpp
  solver/epetra.cpp
  solver/aztecoo.cpp
  solver/amesos.cpp
  solver/mumps.cpp
  solver/superlu.cpp
  solver/petsc.cpp
  solver/umfpack_solver.cpp
  solver/precond_ml.cpp
  solver/precond_ifpack.cpp
  solver/eigensolver.cpp
  solver/eigen.cpp
  compat/fmemopen.cpp
  compat/c99_functions.cpp

  ${TEUCHOS_SOURCES}

  python/python_api.cpp
  python/python_engine.cpp
)

#
# Subdirectories with CMakeLists.txt
#
if(WITH_TESTS)
    enable_testing()
    add_subdirectory(solver)
endif(WITH_TESTS)

if(WITH_PYTHON)
    set(CYTHON_INCLUDE_DIRECTORIES -I${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(hermes_common)
endif(WITH_PYTHON)


# Create file with preprocessor definitions exposing the build settings to the source code.
configure_file(
  config.h.cmake
  config.h
)

#
# Build and install the library.
#

include(BuildAndInstallScripts)
include(PickRealOrCplxLibs)

# This is needed e.g. for eigenvalue problems.
# FIXME: Find which projects need it and enclose their CMake scripts in
# if(WITH_PYTHON) ... endif(WITH_PYTHON). The following should then be not needed,
# as it has been already done in ../CMakeLists.txt
if(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)
    find_package(Python REQUIRED)
    find_package(NumPy REQUIRED)
    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
endif(NOT PYTHON_FOUND OR NOT NUMPY_FOUND)

# Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
# and set its build properties (compile and link flags, installation directories).
macro(BUILD_LIB HERMES_COMMON_LIB BUILD_FLAGS)
  add_library(${HERMES_COMMON_LIB} SHARED ${SRC})
  set_target_properties(${HERMES_COMMON_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})

  PICK_REAL_OR_CPLX_LIBS(${HERMES_COMMON_LIB} ${HERMES_COMMON_LIB})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON_LIB} ${HERMES_COMMON_LIB})

  # Link with all third party libraries corresponding to selected configuration
  # options. Linking of optional libraries for the different Hermeses ({1,2,3}D)
  # is controled in their respective CMakeLists.txt. Note that for unselected
  # TPL's, the respective variables will be empty.
  target_link_libraries(  ${HERMES_COMMON_LIB}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${METIS_LIBRARY}
      ${UMFPACK_LIBRARIES}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_LIBRARIES}
      ${SUPERLU_LIBRARY}
      ${LAPACK_LIBRARIES}
      ${PTHREAD_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
      ${CLAPACK_LIBRARY} ${BLAS_LIBRARY} ${F2C_LIBRARY}
      ${ADDITIONAL_LIBS}
      ${PYTHON_LIBRARY}
  )

  if(NOT MSVC)
      IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
          target_link_libraries(${HERMES_COMMON_LIB} "rt")
      ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      if(HAVE_TEUCHOS_BFD)
          # stacktrace code needs binutils:
          target_link_libraries(${HERMES_COMMON_LIB} bfd iberty)
      endif(HAVE_TEUCHOS_BFD)
  endif(NOT MSVC)
endmacro(BUILD_LIB)

#
# Configure building and installation of all versions of hermes_common needed by
# the selected Hermes configurations.
#

# Build and install real versions:
if(HERMES_COMMON_REAL)
  set(FLAGS "-DHERMES_COMMON_REAL")

  if(MSVC)
    BUILD_LIB(${HERMES_COMMON_LIB_REAL} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_REAL}
                         ${HERMES_COMMON_LIB_REAL_DEBUG}
                         ${HERMES_COMMON_LIB_REAL_RELEASE})
    INSTALL_LIB(${HERMES_COMMON_LIB_REAL})
  else(MSVC)
    if(HERMES_COMMON_DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES_COMMON_LIB_REAL_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_REAL_DEBUG})
    endif(HERMES_COMMON_DEBUG)
    if(HERMES_COMMON_RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES_COMMON_LIB_REAL_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_REAL_RELEASE})
    endif(HERMES_COMMON_RELEASE)
  endif(MSVC)
endif(HERMES_COMMON_REAL)

# Build and install complex versions:
if(HERMES_COMMON_COMPLEX)
  set(FLAGS "-DHERMES_COMMON_COMPLEX")

  if(MSVC)
    BUILD_LIB(${HERMES_COMMON_LIB_CPLX} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_CPLX}
                         ${HERMES_COMMON_LIB_CPLX_DEBUG}
                         ${HERMES_COMMON_LIB_CPLX_RELEASE})
    INSTALL_LIB(${HERMES_COMMON_LIB_CPLX})
  else(MSVC)
    if(HERMES_COMMON_DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES_COMMON_LIB_CPLX_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_CPLX_DEBUG})
    endif(HERMES_COMMON_DEBUG)
    if(HERMES_COMMON_RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES_COMMON_LIB_CPLX_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_CPLX_RELEASE})
    endif(HERMES_COMMON_RELEASE)
  endif(MSVC)
endif(HERMES_COMMON_COMPLEX)


# Install header files.

file(GLOB INC_COMMON    "*.h")
file(GLOB INC_SOLVER    "solver/*.h")
file(GLOB INC_COMPAT    "compat/*.h")
file(GLOB INC_PYTHON    "python/*.h")
file(GLOB INC_TP        "third_party_codes/trilinos-teuchos/*.h*")
install(FILES ${INC_COMMON}    DESTINATION include/hermes/hermes_common)
install(FILES ${INC_SOLVER}    DESTINATION include/hermes/hermes_common/solver)
install(FILES ${INC_COMPAT}    DESTINATION include/hermes/hermes_common/compat)
install(FILES ${INC_PYTHON}    DESTINATION include/hermes/hermes_common/python)
install(FILES ${INC_TP}        DESTINATION include/hermes/hermes_common/third_party_codes/trilinos-teuchos)

