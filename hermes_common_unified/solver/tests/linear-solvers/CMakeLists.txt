project(test-solver)

include(PickRealOrCplxLibs)

if(HERMES_COMMON_REAL)
  add_executable(${PROJECT_NAME} main.cpp)

  if(HERMES_COMMON_DEBUG)
    set(FLAGS "-DHERMES_COMMON_REAL ${DEBUG_FLAGS}")
    set(HERMES_COMMON ${HERMES_COMMON_LIB_REAL_DEBUG})
  else(HERMES_COMMON_DEBUG)
    set(FLAGS "-DHERMES_COMMON_REAL ${RELEASE_FLAGS}")
    set(HERMES_COMMON ${HERMES_COMMON_LIB_REAL_RELEASE})
  endif(HERMES_COMMON_DEBUG)
  
  set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON} ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON} ${TRILINOS_LIBRARIES})
endif(HERMES_COMMON_REAL)

if(HERMES_COMMON_COMPLEX)
  add_executable(${PROJECT_NAME}-cplx main.cpp)
  
  if(HERMES_COMMON_DEBUG)
    set(FLAGS "-DHERMES_COMMON_COMPLEX ${DEBUG_FLAGS}")
    set(HERMES_COMMON ${HERMES_COMMON_LIB_CPLX_DEBUG})
  else(HERMES_COMMON_DEBUG)
    set(FLAGS "-DHERMES_COMMON_COMPLEX ${RELEASE_FLAGS}")
    set(HERMES_COMMON ${HERMES_COMMON_LIB_CPLX_RELEASE})
  endif(HERMES_COMMON_DEBUG)
    
  set_property(TARGET ${PROJECT_NAME}-cplx PROPERTY COMPILE_FLAGS ${FLAGS})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON} ${PROJECT_NAME}-cplx)
  target_link_libraries(${PROJECT_NAME}-cplx ${HERMES_COMMON} ${TRILINOS_LIBRARIES})
endif(HERMES_COMMON_COMPLEX)

# Tests

if(HERMES_COMMON_REAL)
  set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

  if(WITH_PETSC)
    add_test(test-petsc-solver-1 sh -c "${BIN} petsc ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
  endif(WITH_PETSC)

  if(WITH_UMFPACK)
    add_test(test-umfpack-solver-1 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
    add_test(test-umfpack-solver-2 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
    add_test(test-umfpack-solver-3 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

    add_test(test-umfpack-solver-b-1 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
    add_test(test-umfpack-solver-b-2 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
    add_test(test-umfpack-solver-b-3 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
  endif(WITH_UMFPACK)

  if(WITH_TRILINOS)
    if(HAVE_AZTECOO)
      add_test(test-aztecoo-solver-1 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
      add_test(test-aztecoo-solver-2 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
      add_test(test-aztecoo-solver-3 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

      add_test(test-aztecoo-solver-b-1 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
      add_test(test-aztecoo-solver-b-2 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
      add_test(test-aztecoo-solver-b-3 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
    endif(HAVE_AZTECOO)

    if(HAVE_AMESOS)
      add_test(test-amesos-solver-1 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
      add_test(test-amesos-solver-2 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
      add_test(test-amesos-solver-3 sh -c "${BIN} amesos ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

      add_test(test-amesos-solver-b-1 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
      add_test(test-amesos-solver-b-2 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
      add_test(test-amesos-solver-b-3 sh -c "${BIN} amesos-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
    endif(HAVE_AMESOS)
  endif(WITH_TRILINOS)

  if(WITH_MUMPS)
    add_test(test-mumps-solver-1 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
    add_test(test-mumps-solver-2 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
    add_test(test-mumps-solver-3 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")

    add_test(test-mumps-solver-b-1 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
    add_test(test-mumps-solver-b-2 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
    add_test(test-mumps-solver-b-3 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
  endif(WITH_MUMPS)

endif(HERMES_COMMON_REAL)

if(HERMES_COMMON_COMPLEX)
  set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-cplx)

  if(WITH_UMFPACK)
    add_test(test-umfpack-solver-cplx-1 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1") 
    add_test(test-umfpack-solver-cplx-b-1 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
  endif(WITH_UMFPACK)

  if(WITH_TRILINOS)
    if(HAVE_AZTECOO)
      add_test(test-aztecoo-solver-cplx-1 sh -c "${BIN} aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
      add_test(test-aztecoo-solver-cplx-b-1 sh -c "${BIN} aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
    endif(HAVE_AZTECOO)
  endif(WITH_TRILINOS)

  if(WITH_MUMPS)
    add_test(test-mumps-solver-cplx-1 sh -c "${BIN} mumps ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
    add_test(test-mumps-solver-cplx-b-1 sh -c "${BIN} mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-cplx-4 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
  endif(WITH_MUMPS)

endif(HERMES_COMMON_COMPLEX)
