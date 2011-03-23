#include "hermes2d.h"

// This test makes sure that example 32-space-l2 works correctly.

const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Projected function.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform uniform mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an L2 space with default shapeset.
  L2Space space(&mesh, P_INIT);

  // Assemble and solve the finite element problem.
  WeakForm wf_dummy;

  // Initialize the exact and projected solution.
  Solution sln;
  CustomExactSolution sln_exact(&mesh);

  OGProjection::project_global(&space, &sln_exact, &sln, matrix_solver);

  bool success = true;

  if (success == true) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

