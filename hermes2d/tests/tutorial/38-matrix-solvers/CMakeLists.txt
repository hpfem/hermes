project(tutorial-38-matrix-solvers)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

#Possible solvers are: petsc, petsc-block, umfpack, umfpack-block, pardiso, pardiso-block, aztecoo,
#                      aztecoo-block, amesos, amesos-block, mumps,
#                      mumps-block, superlu, superlu-block.
if (WITH_PETSC)
add_test(tutorial-38-matrix-solvers-petsc ${BIN} petsc linsys-1.txt)
add_test(tutorial-38-matrix-solvers-petsc-block ${BIN} petsc-block linsys-1.txt)
endif (WITH_PETSC)

if (WITH_UMFPACK)
add_test(tutorial-38-matrix-solvers-umfpack ${BIN} umfpack linsys-1.txt)
add_test(tutorial-38-matrix-solvers-umfpack-block ${BIN} umfpack-block linsys-1.txt)
endif (WITH_UMFPACK)

if (WITH_PARDISO)
add_test(tutorial-38-matrix-solvers-pardiso ${BIN} pardiso linsys-1.txt)
add_test(tutorial-38-matrix-solvers-pardiso-block ${BIN} pardiso-block linsys-1.txt)
endif (WITH_PARDISO)

if (WITH_TRILINOS)
add_test(tutorial-38-matrix-solvers-aztecoo ${BIN} aztecoo linsys-1.txt)
add_test(tutorial-38-matrix-solvers-aztecoo-block ${BIN} aztecoo-block linsys-1.txt)

add_test(tutorial-38-matrix-solvers-amesos ${BIN} amesos linsys-1.txt)
add_test(tutorial-38-matrix-solvers-amesos-block ${BIN} amesos-block linsys-1.txt)
endif (WITH_TRILINOS)

if (WITH_MUMPS)
add_test(tutorial-38-matrix-solvers-mumps ${BIN} mumps linsys-1.txt)
add_test(tutorial-38-matrix-solvers-mumps-block ${BIN} mumps-block linsys-1.txt)
endif (WITH_MUMPS)

if (WITH_SUPERLU)
add_test(tutorial-38-matrix-solvers-superlu ${BIN} superlu linsys-1.txt)
add_test(tutorial-38-matrix-solvers-superlu-block ${BIN} superlu-block linsys-1.txt)
endif (WITH_SUPERLU)

