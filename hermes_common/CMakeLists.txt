project(hermes_common)
    
  # Check for stacktrace and math library.
  include(CMakeCommon.txt)

  # Additional building steps required by MSVC.
  include(BuildAndInstallScripts)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/solvers)
  
  #
  # Source files for the Hermes_common library.
  #
  set(SRC
    src/mixins.cpp
    src/callstack.cpp
    src/matrix.cpp
    src/cs_matrix.cpp
    src/api.cpp
    src/range.cpp
    src/table.cpp
    src/qsort.cpp
    src/c99_functions.cpp
    src/ord.cpp
    src/hermes_function.cpp
    src/exceptions.cpp
    src/solvers/linear_matrix_solver.cpp
    src/solvers/epetra.cpp
    src/solvers/aztecoo_solver.cpp
    src/solvers/amesos_solver.cpp
    src/solvers/mumps_solver.cpp
    src/solvers/paralution_solver.cpp
    src/solvers/superlu_solver.cpp
    src/solvers/superlu_solver_real.cpp
    src/solvers/superlu_solver_cplx.cpp
    src/solvers/petsc_solver.cpp
    src/solvers/umfpack_solver.cpp
    src/solvers/precond_ml.cpp
    src/solvers/precond_ifpack.cpp
  )

  set(HEADERS
    include/mixins.h
    include/callstack.h
    include/common.h
    include/matrix.h
    include/cs_matrix.h
    include/api.h
    include/array.h
    include/range.h
    include/table.h
    include/qsort.h
    include/c99_functions.h
    include/ord.h
    include/hermes_function.h
    include/exceptions.h
    include/vector.h
    include/solvers/linear_matrix_solver.h
    include/solvers/epetra.h
    include/solvers/aztecoo_solver.h
    include/solvers/amesos_solver.h
    include/solvers/mumps_solver.h
    include/solvers/paralution_solver.h
    include/solvers/superlu_solver.h
    include/solvers/petsc_solver.h
    include/solvers/umfpack_solver.h
    include/solvers/precond_ml.h
    include/solvers/precond_ifpack.h
  )
  
  #
  # Subdirectories with CMakeLists.txt

  # Create file with preprocessor definitions exposing the build settings to the source code.
  configure_file(
    config.h.cmake
    include/config.h
  )

  #
  # Build and install the library.
  #
  include(BuildAndInstallScripts)

  # Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
  # and set its build properties (compile and link flags, installation directories).
  macro(BUILD_COMMON_LIB HERMES_COMMON_LIB BUILD_FLAGS)
    add_library(${HERMES_COMMON_LIB} SHARED ${SRC} ${HEADERS})
    add_definitions(-DHERMES_COMMON)
    set_target_properties(${HERMES_COMMON_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})

    # Link with all third party libraries corresponding to selected configuration
    # options. Linking of optional libraries for the different Hermeses ({1, 2, 3}D)
    # is controled in their respective CMakeLists.txt. Note that for unselected
    # TPL's, the respective variables will be empty.

    if(NOT(${WITH_TC_MALLOC}))
      set(TCMALLOC_LIBRARY "")
    endif(NOT(${WITH_TC_MALLOC}))
    
    target_link_libraries(  ${HERMES_COMMON_LIB}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${METIS_LIBRARY}
      ${UMFPACK_LIBRARIES} ${PARALUTION_LIBRARY}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_LIBRARIES}
      ${SUPERLU_LIBRARY}
      ${LAPACK_LIBRARIES}
      ${STACK_WALKER_LIBRARY}
      ${PTHREAD_LIBRARY} ${TCMALLOC_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
      ${CLAPACK_LIBRARY} ${BLAS_LIBRARY} ${F2C_LIBRARY}
      ${WINBLAS_LIBRARY} 
      ${ADDITIONAL_LIBS}
    )

    if(NOT MSVC)
      IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        target_link_libraries(${HERMES_COMMON_LIB} "rt")
      ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    endif(NOT MSVC)
  endmacro(BUILD_COMMON_LIB)

  #
  # Configure building and installation of all versions of hermes_common needed by
  # the selected Hermes configurations.
  #
  if(MSVC)
    if(${HERMES_COMMON_DEBUG})
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_DEBUG} ${DEBUG_FLAGS})
      ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_DEBUG}
                         ${HERMES_COMMON_LIB_DEBUG}
                         ${HERMES_COMMON_LIB_RELEASE})
      INSTALL_LIB(${HERMES_COMMON_LIB_DEBUG})
      set(HERMES_COMMON ${HERMES_COMMON_LIB_DEBUG})
    endif(${HERMES_COMMON_DEBUG})
    
    if(${HERMES_COMMON_RELEASE})
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_RELEASE} ${RELEASE_FLAGS})
      ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_RELEASE}
                         ${HERMES_COMMON_LIB_RELEASE}
                         ${HERMES_COMMON_LIB_RELEASE})
      INSTALL_LIB(${HERMES_COMMON_LIB_RELEASE})
      set(HERMES_COMMON ${HERMES_COMMON_LIB_RELEASE})
    endif(${HERMES_COMMON_RELEASE})
  else(MSVC)
    if(HERMES_COMMON_DEBUG)
      set(BUILD_FLAGS "${DEBUG_FLAGS}")
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_DEBUG})
    endif(HERMES_COMMON_DEBUG)
    
    if(HERMES_COMMON_RELEASE)
      set(BUILD_FLAGS "${RELEASE_FLAGS}")
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_RELEASE})
    endif(HERMES_COMMON_RELEASE)
  endif(MSVC)

  # Install header files.

  file(GLOB INC_COMMON    "include/*.h")
  file(GLOB INC_SOLVER    "include/solvers/*.h")
  install(FILES ${INC_COMMON}    DESTINATION ${TARGET_ROOT}/include/hermes_common)
  install(FILES ${INC_SOLVER}    DESTINATION ${TARGET_ROOT}/include/hermes_common/solvers)
