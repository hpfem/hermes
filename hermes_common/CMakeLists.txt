cmake_minimum_required(VERSION 2.6)

project(hermes_common)

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
add_subdirectory(solver)

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

# Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
# and set its build properties (compile and link flags, installation directories).
macro(BUILD_LIB HERMES_COMMON_LIB BUILD_FLAGS)
  add_library(${HERMES_COMMON_LIB} SHARED ${SRC})
  set_target_properties(${HERMES_COMMON_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})

  #This is needed e.g. for eigenvalue problems.
  if(NOT Python_FOUND OR NOT NumPy_FOUND)
    find_package(Python REQUIRED)
    find_package(NumPy REQUIRED)
    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
    
    # This path will be used relatively to CMAKE_INSTALL_PREFIX.
    set(PYTHON_INSTALL_PATH lib/python/site-packages)   
  endif(NOT Python_FOUND OR NOT NumPy_FOUND)    

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
install(FILES ${INC_COMMON}    DESTINATION include/hermes/hermes_common)
install(FILES ${INC_SOLVER}    DESTINATION include/hermes/hermes_common/solver)
install(FILES ${INC_COMPAT}    DESTINATION include/hermes/hermes_common/compat)


