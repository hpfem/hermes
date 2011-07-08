project(test-solvers-complex)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON_LIB} ${TRILINOS_LIBRARIES})

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

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