Using Matrix Solvers (40-matrix-solvers)
-------------------------

**Git reference:** Tutorial example `40-matrix-solvers <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P10-miscellaneous/40-matrix-solvers>`_. 

This example shows how to solve a matrix problem Ax = b in Hermes using various solvers
that Hermes supports. The matrix and the right-hand side vector are first read from a text 
file into coordinate matrix structures, then passed into a SparseMatrix and Vector structures, 
and finally solved using a LinearSolver of the user's choice.

To use a solver, Hermes2D must be compiled with this solver enabled. This is done via
compilation flags WITH_PETSC, WITH_UMFPACK, WITH_TRILINOS (for AztecOO and Amesos),
WITH_MUMPS and WITH_SUPERLU that are set to default values in CMakeLists.txt. To change these settings, 
use the CMake.vars file. For example, a line "set(WITH_MUMPS YES)" will override the 
default value for the Mumps solver.

**Usage:**
::

    matrix-solvers <matrix_solver> <input_text_file_with_matrix_and_vector> 

**Available solvers:** petsc, petsc-block, umfpack, umfpack-block,  aztecoo, 
aztecoo-block, amesos, amesos-block, mumps, mumps-block, superlu, superlu-block.

**Available input files:** linsys-1.txt, linsys-2.txt, linsys-3.txt,

but the user can of course provide his/her own. In the input file, the first number is the 
matrix/vector size, followed by the matrix in the coordinate format and the right-hand side 
vector in the coordinate format. For example the file linsys-3.txt contains::

    5
    0 0 2
    0 1 3
    1 0 3
    1 2 4
    1 4 6
    2 1 -1
    2 2 -3
    2 3 2
    3 2 1
    4 1 4
    4 2 2
    4 4 1

    0 8
    1 45
    2 -3
    3 3
    4 19

Output::

    I Setting matrix element (0, 0, 2).
    I Setting matrix element (0, 1, 3).
    I Setting matrix element (1, 0, 3).
    I Setting matrix element (1, 2, 4).
    I Setting matrix element (1, 4, 6).
    I Setting matrix element (2, 1, -1).
    I Setting matrix element (2, 2, -3).
    I Setting matrix element (2, 3, 2).
    I Setting matrix element (3, 2, 1).
    I Setting matrix element (4, 1, 4).
    I Setting matrix element (4, 2, 2).
    I Setting matrix element (4, 4, 1).
    I Setting vector entry (0, 8).
    I Setting vector entry (1, 45).
    I Setting vector entry (2, -3).
    I Setting vector entry (3, 3).
    I Setting vector entry (4, 19).
    I Matrix solve successful.
    Solution vector: 1 2 3 4 5 
