project(hermes_common)
    
  # Check for stacktrace and math library.
  include(CheckFunctionExists)
  CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
  CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)

  # Additional building steps required by MSVC.
  include(BuildAndInstallScripts)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/algebra)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/util)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/solvers)
  
  #
  # Source files for the Hermes_common library.
  #
  set(SRC
    src/mixins.cpp
    src/ord.cpp
    src/hermes_function.cpp
    src/exceptions.cpp
    src/api.cpp
    src/algebra/algebra_utilities.cpp
    src/algebra/matrix.cpp
    src/algebra/vector.cpp
    src/algebra/algebra_mixins.cpp
    src/algebra/dense_matrix_operations.cpp
    src/algebra/cs_matrix.cpp
    src/util/callstack.cpp
    src/util/qsort.cpp
    src/util/c99_functions.cpp
    src/data_structures/range.cpp
    src/data_structures/table.cpp
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
    include/common.h
    include/api.h
    include/ord.h
    include/hermes_function.h
    include/exceptions.h
    include/util/c99_functions.h
    include/util/compat.h
    include/util/callstack.h
    include/util/qsort.h
    include/algebra/algebra_utilities.h
    include/algebra/matrix.h
    include/algebra/vector.h
    include/algebra/cs_matrix.h
    include/algebra/algebra_mixins.h
    include/algebra/dense_matrix_operations.h
    include/data_structures/array.h
    include/data_structures/range.h
    include/data_structures/table.h
    include/data_structures/hermes_vector.h
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
    include/solvers/precond.h
  )
  
  SOURCE_GROUP(
    "Source Files\\Matrix Solvers" FILES 
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
  
  SOURCE_GROUP(
    "Source Files\\Algebra" FILES 
    src/algebra/algebra_utilities.cpp
    src/algebra/matrix.cpp
    src/algebra/vector.cpp
    src/algebra/algebra_mixins.cpp
    src/algebra/dense_matrix_operations.cpp
    src/algebra/cs_matrix.cpp
  )
  
  SOURCE_GROUP(
    "Source Files\\Data Structures" FILES 
    src/data_structures/range.cpp
    src/data_structures/table.cpp
  )
  
  SOURCE_GROUP(
    "Source Files\\Utilities" FILES 
    src/util/callstack.cpp
    src/util/qsort.cpp
    src/util/c99_functions.cpp
  )
  
  SOURCE_GROUP(
    "Header Files\\Matrix Solvers" FILES 
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
    include/solvers/precond.h
  )
  
  SOURCE_GROUP(
    "Header Files\\Algebra" FILES 
    include/algebra/algebra_utilities.h
    include/algebra/matrix.h
    include/algebra/vector.h
    include/algebra/cs_matrix.h
    include/algebra/algebra_mixins.h
    include/algebra/dense_matrix_operations.h
  )
  
  SOURCE_GROUP(
    "Header Files\\Data Structures" FILES 
    include/data_structures/array.h
    include/data_structures/range.h
    include/data_structures/table.h
    include/data_structures/hermes_vector.h
  )
  
  SOURCE_GROUP(
    "Header Files\\Utilities" FILES 
    include/util/c99_functions.h
    include/util/compat.h
    include/util/callstack.h
    include/util/qsort.h
  )
  
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
    if(${HERMES_STATIC_LIBS})
      add_library(${HERMES_COMMON_LIB} ${SRC} ${HEADERS})
    else()
      add_library(${HERMES_COMMON_LIB} SHARED ${SRC} ${HEADERS})
    endif()
    
    add_definitions(-DHERMES_COMMON)
    if(NOT MSVC)
      set_target_properties(${HERMES_COMMON_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
    endif()
    
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
      ${PTHREAD_LIBRARY}
      ${TCMALLOC_LIBRARY}
      ${MPI_LIBRARIES}
      ${SCALAPACK_LIBRARIES}
      ${CLAPACK_LIBRARY}
      ${BLAS_LIBRARY}
      ${F2C_LIBRARY}
      ${BSON_LIBRARY}
      ${MATIO_LIBRARY}
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
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_DEBUG} ${HERMES_DEBUG_FLAGS})
      ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_DEBUG}
                         ${HERMES_COMMON_LIB_DEBUG}
                         ${HERMES_COMMON_LIB_RELEASE})
      INSTALL_LIB(${HERMES_COMMON_LIB_DEBUG})
      set(HERMES_COMMON ${HERMES_COMMON_LIB_DEBUG})
    endif(${HERMES_COMMON_DEBUG})
    
    if(${HERMES_COMMON_RELEASE})
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_RELEASE} ${HERMES_RELEASE_FLAGS})
      ADD_MSVC_BUILD_FLAGS(${HERMES_COMMON_LIB_RELEASE}
                         ${HERMES_COMMON_LIB_RELEASE}
                         ${HERMES_COMMON_LIB_RELEASE})
      INSTALL_LIB(${HERMES_COMMON_LIB_RELEASE})
      set(HERMES_COMMON ${HERMES_COMMON_LIB_RELEASE})
    endif(${HERMES_COMMON_RELEASE})
  else(MSVC)
    if(HERMES_COMMON_DEBUG)
      set(BUILD_FLAGS "${HERMES_DEBUG_FLAGS}")
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_DEBUG})
    endif(HERMES_COMMON_DEBUG)
    
    if(HERMES_COMMON_RELEASE)
      set(BUILD_FLAGS "${HERMES_RELEASE_FLAGS}")
      BUILD_COMMON_LIB(${HERMES_COMMON_LIB_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES_COMMON_LIB_RELEASE})
    endif(HERMES_COMMON_RELEASE)
  endif(MSVC)

  # Install header files.

  file(GLOB INC_COMMON    "include/*.h")
  file(GLOB INC_SOLVER    "include/solvers/*.h")
  file(GLOB INC_STRUCTURES "include/data_structures/*.h")
  file(GLOB INC_ALGEBRA    "include/algebra/*.h")
  file(GLOB INC_UTIL    "include/util/*.h")
  install(FILES ${INC_COMMON}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hermes_common)
  install(FILES ${INC_SOLVER}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hermes_common/solvers)
  install(FILES ${INC_STRUCTURES}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hermes_common/data_structures)
  install(FILES ${INC_ALGEBRA}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hermes_common/algebra)
  install(FILES ${INC_UTIL}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hermes_common/util)
