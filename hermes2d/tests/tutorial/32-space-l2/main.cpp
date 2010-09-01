#include "hermes2d.h"

// This test makes sure that example 32-space-l2 works correctly.

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

  // Assemble and solve the finite element problem.
  WeakForm wf_dummy;
  LinearProblem ls(&wf_dummy, &space);
  Solution* sln_tmp = new Solution(&mesh, F);
  Solution sln;
  project_global(&space, H2D_L2_NORM, sln_tmp, &sln);
  delete sln_tmp;

  // Visualize the solution.
  ScalarView view1("Projection", 610, 0, 600, 500);

  // It will "Exception: SegFault" if we do not use View::wait() or View::close(). 
  view1.close();

  bool success = true;

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success == true) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

