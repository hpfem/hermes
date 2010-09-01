#include "hermes2d.h"

// This example shows how to use the L2 finite element space and L2 shapeset.
// As a sample problem, a continuous function x^3 + y^3 is projected onto the
// L2 finite element space in the L2 norm. When zero-order is used, the result
// is a piecewice constant function. The class BaseView will show you the basis
// functions.
//
// The following parameters can be changed:

const int INIT_REF_NUM = 1;    // Number of initial uniform mesh refinements.
const int P_INIT = 3;          // Polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

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

  // Create an L2 space with default shapeset.
  L2Space space(&mesh, P_INIT);

  // View basis functions.
  BaseView bview("BaseView", 0, 0, 600, 500);
  bview.show(&space);
  View::wait(H2DV_WAIT_KEYPRESS);

  // Assemble and solve the finite element problem.
  WeakForm wf_dummy;
  LinearProblem ls(&wf_dummy, &space);
  Solution* sln_tmp = new Solution(&mesh, F);
  Solution sln;
  project_global(&space, H2D_L2_NORM, sln_tmp, &sln);
  delete sln_tmp;

  // Visualize the solution.
  ScalarView view1("Projection", 610, 0, 600, 500);
  view1.show(&sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

