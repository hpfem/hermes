project(test-solvers-complex)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON_LIB} ${TRILINOS_LIBRARIES})

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

if(WITH_PETSC)
  add_test(test-petsc-solver-cplx-1 ${BIN} petsc)
  add_test(test-petsc-solver-cplx-b-1 ${BIN} petsc-block)
endif(WITH_PETSC)

if(WITH_UMFPACK)
  add_test(test-umfpack-solver-cplx-1 ${BIN} umfpack) 
  add_test(test-umfpack-solver-cplx-b-1 ${BIN} umfpack-block)
endif(WITH_UMFPACK)

if(WITH_TRILINOS)
  if(HAVE_AMESOS)
    add_test(test-amesos-solver-cplx-1 ${BIN} amesos)
    add_test(test-amesos-solver-cplx-b-1 ${BIN} amesos-block)
  endif(HAVE_AMESOS)
endif(WITH_TRILINOS)

if(WITH_TRILINOS)
  if(HAVE_AZTECOO)
    add_test(test-aztecoo-solver-cplx-1 ${BIN} aztecoo)
    add_test(test-aztecoo-solver-cplx-b-1 ${BIN} aztecoo-block)
  endif(HAVE_AZTECOO)
endif(WITH_TRILINOS)

if(WITH_MUMPS)
  add_test(test-mumps-solver-cplx-1 ${BIN} mumps)
  add_test(test-mumps-solver-cplx-b-1 ${BIN} mumps-block)
endif(WITH_MUMPS)
