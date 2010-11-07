project(solver)

include (${hermes3d_SOURCE_DIR}/CMake.common)
	
if(H3D_REAL)
	add_executable(${PROJECT_NAME}		main.cpp)
  set_common_target_properties(${PROJECT_NAME} ${HERMES_REAL_BIN})
endif(H3D_REAL)
if(H3D_COMPLEX)
	add_executable(${PROJECT_NAME}-cplx		main.cpp)
  set_common_target_properties(${PROJECT_NAME}-cplx ${HERMES_CPLX_BIN})
endif(H3D_COMPLEX)

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
if(H3D_REAL)
if(WITH_PETSC)
	add_test(petsc-solver-1 "${BIN}" petsc "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
#	add_test(petsc-solver-2 "${BIN}" petsc "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
endif(WITH_PETSC)

if(WITH_UMFPACK)
	add_test(umfpack-solver-1 "${BIN}" umfpack "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(umfpack-solver-2 "${BIN}" umfpack "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(umfpack-solver-3 "${BIN}" umfpack "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	add_test(umfpack-solver-4 "${BIN}" umfpack "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
	
	add_test(umfpack-solver-b-1 "${BIN}" umfpack-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(umfpack-solver-b-2 "${BIN}" umfpack-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(umfpack-solver-b-3 "${BIN}" umfpack-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	add_test(umfpack-solver-b-4 "${BIN}" umfpack-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
endif(WITH_UMFPACK)

if(WITH_PARDISO)
	add_test(pardiso-solver-1 "${BIN}" pardiso "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(pardiso-solver-2 "${BIN}" pardiso "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(pardiso-solver-3 "${BIN}" pardiso "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	# NOTE: PARDISO is not able to recognize singular system	
#	add_test(pardiso-solver-4 "${BIN}" pardiso "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
	
	add_test(pardiso-solver-b-1 "${BIN}" pardiso-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(pardiso-solver-b-2 "${BIN}" pardiso-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(pardiso-solver-b-3 "${BIN}" pardiso-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	# NOTE: PARDISO is not able to recognize singular system	
#	add_test(pardiso-solver-b-4 "${BIN}" pardiso-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
endif(WITH_PARDISO)

if(WITH_TRILINOS)
  if(HAVE_AZTECOO)
	  add_test(aztecoo-solver-1 "${BIN}" aztecoo "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	  add_test(aztecoo-solver-2 "${BIN}" aztecoo "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	  add_test(aztecoo-solver-3 "${BIN}" aztecoo "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")

	  add_test(aztecoo-solver-b-1 "${BIN}" aztecoo-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	  add_test(aztecoo-solver-b-2 "${BIN}" aztecoo-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	  add_test(aztecoo-solver-b-3 "${BIN}" aztecoo-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
  endif(HAVE_AZTECOO)
  
  if(HAVE_AMESOS)
	  add_test(amesos-solver-1 "${BIN}" amesos "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	  add_test(amesos-solver-2 "${BIN}" amesos "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	  add_test(amesos-solver-3 "${BIN}" amesos "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")

	  add_test(amesos-solver-b-1 "${BIN}" amesos-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	  add_test(amesos-solver-b-2 "${BIN}" amesos-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	  add_test(amesos-solver-b-3 "${BIN}" amesos-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	endif(HAVE_AMESOS)
endif(WITH_TRILINOS)

if(WITH_MUMPS)
	add_test(mumps-solver-1 "${BIN}" mumps "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(mumps-solver-2 "${BIN}" mumps "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(mumps-solver-3 "${BIN}" mumps "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	add_test(mumps-solver-4 "${BIN}" mumps "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
	
	add_test(mumps-solver-b-1 "${BIN}" mumps-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1")
	add_test(mumps-solver-b-2 "${BIN}" mumps-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2")
	add_test(mumps-solver-b-3 "${BIN}" mumps-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3")
	add_test(mumps-solver-b-4 "${BIN}" mumps-block "${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular")
endif(WITH_MUMPS)

endif(H3D_REAL)

if(H3D_COMPLEX)
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-cplx)

if(WITH_UMFPACK)
	add_test(umfpack-solver-cplx-1 "${BIN}" umfpack ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	add_test(umfpack-solver-cplx-b-1 "${BIN}" umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
endif(WITH_UMFPACK)

if(WITH_TRILINOS)
  if(HAVE_AZTECOO)
	  add_test(aztecoo-solver-cplx-1 "${BIN}" aztecoo ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	  add_test(aztecoo-solver-cplx-b-1 "${BIN}" aztecoo-block ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	endif(HAVE_AZTECOO)
endif(WITH_TRILINOS)

if(WITH_MUMPS)
	add_test(mumps-solver-cplx-1 "${BIN}" mumps ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
	add_test(mumps-solver-cplx-b-1 "${BIN}" mumps-block ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-cplx-1")
endif(WITH_MUMPS)

endif(H3D_COMPLEX)
