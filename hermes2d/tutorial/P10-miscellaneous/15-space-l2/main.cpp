#include "hermes2d.h"

// This example shows how to use the L2 finite element space and L2 shapeset.
// As a sample problem, a continuous function x^3 + y^3 is projected onto the
// L2 finite element space in the L2 norm. When zero-order is used, the result
// is a piecewice constant function. The class BaseView will show you the basis
// functions.
//
// The following parameters can be changed:

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
  for (int i = 0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers (default is Neumann boundary).
  BCTypes bc_types;

  // Enter Dirichlet boundary values (not relevant here).
  BCValues bc_values;

  // Create an L2 space with default shapeset.
  L2Space space(&mesh, &bc_types, &bc_values, P_INIT);

  // View basis functions.
  BaseView bview("BaseView", new WinGeom(0, 0, 600, 500));
  bview.show(&space);
  // View::wait(H2DV_WAIT_KEYPRESS);

  // Assemble and solve the finite element problem.
  WeakForm wf_dummy;

  // Initialize the exact and projected solution.
  Solution sln;
  Solution sln_exact(&mesh, F);

  OGProjection::project_global(&space, &sln_exact, &sln, matrix_solver);

  // Visualize the solution.
  ScalarView view1("Projection", new WinGeom(610, 0, 600, 500));
  view1.show(&sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

