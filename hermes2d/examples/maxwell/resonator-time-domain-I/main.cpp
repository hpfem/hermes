#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example shows how to discretize the first-order time-domain Maxwell's equations 
// with vector-valued E (an Hcurl function) and scalar B (an L2 function). Time integration 
// is done using implicit Euler.
//
// PDE: (1. / SPEED_OF_LIGHT**2) \partial E / \partial t - curl B = 0,
//      \partial B / \partial t + curl E = 0.
//
// Note: curl E = \partial E_2 / \partial x - \partial E_1 / \partial y
//       curl B = (\partial B / \partial y, - \partial B / \partial x)
//       
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: Rectangular domain (-a, a) x (-b, b)... See mesh file domain.mesh.
//
// BC: perfect conductor for E on the entire boundary.
//
// The following parameters can be changed:

const int P_INIT_E = 1;                           // Initial polynomial degree for E.
const int P_INIT_B = 1;                           // Initial polynomial degree for B.
const int INIT_REF_NUM = 5;                       // Number of initial uniform mesh refinements.
const double time_step = 0.005;                   // Time step.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double SPEED_OF_LIGHT = 1.0;                // Speed of light.
const double T_FINAL = 10000;                     // Length of time interval.
double current_time = 0;

// Initial conditions for E and B (Version A).
scalar2 E_init_cond(double x, double y, scalar2& dx, scalar2& dy) {
  dx[0] = sin(x) * sin(y);
  dx[1] = cos(x) * cos(y);
  dy[0] = -cos(x) * cos(y);
  dy[1] = -sin(x) * sin(y);
  return scalar2(-cos(x) * sin(y), sin(x) * cos(y));
}
scalar B_init_cond(double x, double y, scalar& dx, scalar& dy) {
  dx = 0.0;
  dy = 0.0;
  return 0.0;
}

/* 
// Initial conditions for E and B (Version B).
scalar2 E_init_cond(double x, double y, scalar2& dx, scalar2& dy) {
  dx[0] = 0.0;
  dx[1] = 0.0;
  dy[0] = 0.0;
  dy[1] = 0.0;
  return scalar2(0.0, 0.0);
}
scalar B_init_cond(double x, double y, scalar& dx, scalar& dy) {
  dx = - 0.5 * M_PI * sin(0.5 * M_PI * x) * cos(0.5 * M_PI * y);
  dy = - 0.5 * M_PI * cos(0.5 * M_PI * x) * sin(0.5 * M_PI * y);
  return cos(0.5 * M_PI * x) * cos(0.5 * M_PI * y);
}
*/

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

  // Initialize boundary conditions.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create an Hcurl space for E and L2 space for B.
  HcurlSpace E_space(&mesh, &bc_types, &bc_values, P_INIT_E);
  L2Space B_space(&mesh, P_INIT_B);
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&E_space, &B_space));
  info("ndof = %d", ndof);

  // Initialize solutions for E and B.
  Solution E_sln(&mesh, E_init_cond);
  Solution B_sln(&mesh, B_init_cond);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, &E_sln);
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, &B_sln);

  // Initialize views.
  ScalarView E1_view("E1", new WinGeom(0, 0, 520, 400));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("E2", new WinGeom(525, 0, 520, 400));
  E2_view.fix_scale_width(50);
  ScalarView B_view("B", new WinGeom(1050, 0, 520, 400));
  B_view.fix_scale_width(50);
  
  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&E_space, &B_space), is_linear);

  // Initialize matrix solver.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Time stepping:
  int ts = 1; bool rhs_only = false;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
    rhs_only = true;

    // Solve the linear system of the reference problem. 
    // If successful, obtain the solutions.
    info("Solving.");
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), 
                                  Hermes::vector<Space*>(&E_space, &B_space), 
				  Hermes::vector<Solution*>(&E_sln, &B_sln));
    else error ("Matrix solver failed.\n");

    // Show the solution.
    E1_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    E2_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);
    B_view.show(&B_sln);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;

  } while (current_time < T_FINAL);

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

