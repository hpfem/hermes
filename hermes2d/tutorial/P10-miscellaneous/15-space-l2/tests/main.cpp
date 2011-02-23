#include "hermes2d.h"

// This test makes sure that example 32-space-l2 works correctly.

const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Projected function.
scalar F(double x, double y, double& dx, double& dy)
{
  return - pow(x, 4) * pow(y, 5); 
  dx = 0; // not needed for L2-projection
  dy = 0; // not needed for L2-projection
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform uniform mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;

  // Enter Dirichlet boundary values.
  BCValues bc_values;

  // Create an L2 space with default shapeset.
  L2Space space(&mesh, &bc_types, &bc_values, P_INIT);


  // Assemble and solve the finite element problem.
  WeakForm wf_dummy;

  // Initialize the exact and projected solution.
  Solution sln;
  Solution sln_exact(&mesh, F);

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

