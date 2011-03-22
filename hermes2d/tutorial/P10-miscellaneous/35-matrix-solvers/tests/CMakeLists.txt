project(test-tutorial-P10-35-matrix-solvers)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

#Possible solvers are: petsc, petsc-block, umfpack, umfpack-block, aztecoo, aztecoo-block, 
#                      amesos, amesos-block, mumps, mumps-block, superlu, superlu-block.
if (WITH_PETSC)
add_test(test-tutorial-P10-35-matrix-solvers-petsc ${BIN} petsc ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-petsc-block ${BIN} petsc-block ../in/linsys-1.txt)
endif (WITH_PETSC)

if (WITH_UMFPACK)
add_test(test-tutorial-P10-35-matrix-solvers-umfpack ${BIN} umfpack ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-umfpack-block ${BIN} umfpack-block ../in/linsys-1.txt)
endif (WITH_UMFPACK)

if (WITH_TRILINOS)
add_test(test-tutorial-P10-35-matrix-solvers-aztecoo ${BIN} aztecoo ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-aztecoo-block ${BIN} aztecoo-block ../in/linsys-1.txt)

add_test(test-tutorial-P10-35-matrix-solvers-amesos ${BIN} amesos ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-amesos-block ${BIN} amesos-block ../in/linsys-1.txt)
endif (WITH_TRILINOS)

if (WITH_MUMPS)
add_test(test-tutorial-P10-35-matrix-solvers-mumps ${BIN} mumps ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-mumps-block ${BIN} mumps-block ../in/linsys-1.txt)
endif (WITH_MUMPS)

if (WITH_SUPERLU)
add_test(test-tutorial-P10-35-matrix-solvers-superlu ${BIN} superlu ../in/linsys-1.txt)
add_test(test-tutorial-P10-35-matrix-solvers-superlu-block ${BIN} superlu-block ../in/linsys-1.txt)
endif (WITH_SUPERLU)

