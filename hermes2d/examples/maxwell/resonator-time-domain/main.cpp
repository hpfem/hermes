#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example shows how to discretize time-domain Maxwell system with vector-valued
// E (an Hcurl function) and scalar B (an L2 function). Time integration is done using 
// implicit Euler. The speed of light is assumed to be 1.0.
//
// PDE: \partial E / \partial t - curl B = 0,
//      \partial B / \partial t + curl E = 0.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: Rectangular domain (-a, a) x (-b, b)... See mesh file domain.mesh.
//
// BC: perfect conductor for E on the entire boundary.
//
// The following parameters can be changed:

const double TAU = 0.1;                           // Time step.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree. NOTE: The meaning is different from
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Initial condition for E.
scalar2 E_init_cond(double x, double y, scalar2& dx, scalar2& dy) {
  dx[0] = 0.0;
  dx[1] = -2*x;
  dy[0] = -2*y;
  dy[1] = 0.0;
  return scalar2(1 - y*y, 1 - x*x);
}

//  Boundary markers.
const int BDY_DIRICHLET = 1;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create an Hcurl space for E and L2 space for B.
  HcurlSpace E_space(&mesh, &bc_types, &bc_values, P_INIT);
  L2Space B_space(&mesh, P_INIT);

  // Initialize previous time level solutions;
  Solution E_prev(&mesh, E_init_cond);
  Solution B_prev(&mesh, 0.0);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, &E_prev);
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, &B_prev);

  // Initialize views.
  VectorView E_view("Electric field", new WinGeom(0, 0, 580, 400));
  ScalarView B_view("Magnetic field", new WinGeom(590, 0, 580, 400));
  
  // Initialize matrix solver.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Assemble the linear system.
  info("Assembling.");
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(&wf, Hermes::vector<Space*>(&E_space, &B_space), is_linear);
  dp->assemble(matrix, rhs);
    
  // Solve the linear system of the reference problem. 
  // If successful, obtain the solution.
  Solution E_new(&mesh);
  Solution B_new(&mesh);
  if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), 
                                Hermes::vector<Space*>(&E_space, &B_space), 
				Hermes::vector<Solution*>(&E_new, &B_new));
  else error ("Matrix solver failed.\n");

  // Show real part of the solution.
  E_view.show(&E_new);
  B_view.show(&B_new);

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

