#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to solve a matrix problem Ax = b in Hermes. 
// The matrix and the right-hand side vector are first read from a file 
// in the coordinate format, then passed into a SparseMatrix and Vector,
// and then solved using a LinearSolver.
//
// Possible solvers: petsc, petsc-block, umfpack, umfpack-block, aztecoo, aztecoo-block, 
//                   amesos, amesos-block, mumps, mumps-block, superlu, superlu-block.
//
// Sample usage: "matrix-solvers umfpack linsys-N.txt" where N = 1, 2, 3.

// Include helpers.
#include "../helpers.cpp"

// Main function
int main(int argc, char *argv[]) {

  // PETSc initialization, if Hermes is compiled with PETSc support.
# ifdef WITH_PETSC
  // Do NOT forget to call this when using PETSc solver.
  PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
  // Disable PETSc error handler.
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);
# endif

  // Check number of command-line parameters.
  if (argc < 3) {
    warn("Possible solvers are: petsc, petsc-block, umfpack, umfpack-block, \
aztecoo, aztecoo-block, amesos, amesos-block, mumps, mumps-block, superlu, superlu-block.");
    error("Not enough parameters: Provide a solver and an input file with a matrix and vector.");
  }

  // Prepare to read from file.
  int n;                               // Matrix size.
  int nnz;                             // Number on non-zero entries
  bool cplx_2_real;                    // Decides do we turn complex matrix to real
  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, scalar> ar_rhs; 

  if (argc == 4 && strcasecmp(argv[3],"complex-matrix-to-real") == 0)
     cplx_2_real = true;
  else
     cplx_2_real = false;

  // Read matrix and rhs from file.
  info("Reading matrix from file %s.", argv[2]);
  if (!read_matrix_and_rhs(argv[2], n, nnz, ar_mat, ar_rhs, cplx_2_real))
    error("Failed to read the matrix and rhs.");

  if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    PetscMatrix mat;
    PetscVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("PETSc solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with PETSc support.");
#endif
  }
  else if (strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    PetscMatrix mat;
    PetscVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("PETSc solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with PETSc support.");
#endif
  }
  else if (strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix mat;
    UMFPackVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("UMFpack solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with UMFpack support.");
#endif
  }
  else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix mat;
    UMFPackVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("UMFpack solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with UMFpack support.");
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo") == 0) {
#if defined WITH_TRILINOS && defined HAVE_AZTECOO
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("AztecOO solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with Trilinos support.");
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo-block") == 0) {
#if defined WITH_TRILINOS && defined HAVE_AZTECOO
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("AztecOO solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with Trilinos support.");
#endif
  }
  else if (strcasecmp(argv[1], "amesos") == 0) {
#if defined WITH_TRILINOS && defined HAVE_AMESOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver::is_available("Klu")) {
      AmesosSolver solver("Klu", &mat, &rhs);
      if (solver.solve()) {
        scalar *sln = solver.get_solution();
        info("Matrix solve successful.");
        printf("Solution vector: ");
        if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        printf("\n");
      }
      else {
        info("Amesos solver failed.");
        return ERR_FAILURE;
      }
    }
    else {
      info("Amesos solver failed - Klu not available.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with Trilinos support.");
#endif
  }
  else if (strcasecmp(argv[1], "amesos-block") == 0) {
#if defined WITH_TRILINOS && defined HAVE_AMESOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver::is_available("Klu")) {
      AmesosSolver solver("Klu", &mat, &rhs);
      if (solver.solve()) {
        scalar *sln = solver.get_solution();
        info("Matrix solve successful.");
        printf("Solution vector: ");
        if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
        else {printf("Failed. \n");  return ERR_FAILURE;} 
        printf("\n");
      }
      else {
        info("Amesos solver failed.");
        return ERR_FAILURE;
      }
    }
    else {
      info("Amesos solver failed - Klu not available.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with Amesos support.");
#endif
  }
  else if (strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix mat;
    MumpsVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    } 
    else {
      info("Mumps solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with MUMPS support.");
#endif
  }
  else if (strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix mat;
    MumpsVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("Mumps solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with MUMPS support.");
#endif
  }  
  else if (strcasecmp(argv[1], "superlu") == 0) {
#ifdef WITH_SUPERLU
    SuperLUMatrix mat;
    SuperLUVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    SuperLUSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("SuperLU solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with SuperLU support.");
#endif
  }
  else if (strcasecmp(argv[1], "superlu-block") == 0) {
#ifdef WITH_SUPERLU
    SuperLUMatrix mat;
    SuperLUVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    SuperLUSolver solver(&mat, &rhs);
    if (solver.solve()) {
      scalar *sln = solver.get_solution();
      info("Matrix solve successful.");
      printf("Solution vector: ");
      if ((sln[0] - 4) < 1e-6) printf("%g ", sln[0]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[1] - 2) < 1e-6) printf("%g ", sln[1]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      if ((sln[2] - 3) < 1e-6) printf("%g ", sln[2]);
      else {printf("Failed. \n");  return ERR_FAILURE;} 
      printf("\n");
    }
    else {
      info("SuperLU solver failed.");
      return ERR_FAILURE;
    }
#else
    error("Hermes was not built with SuperLU support.");
#endif
  }  
  else {
    warn("Possible solvers are: petsc, petsc-block, umfpack, umfpack-block, \
aztecoo, aztecoo-block, amesos, amesos-block, mumps, mumps-block, superlu, superlu-block.");
    error("Unknown matrix solver.");
  }

#ifdef WITH_PETSC
  // Do NOT forget to call this when using PETSc solver.
  PetscFinalize();
#endif


  return ERR_SUCCESS;
}
