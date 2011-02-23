#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example uses the Newton's method to solve a nonlinear complex-valued
//  time-dependent PDE (the Gross-Pitaevski equation describing the behavior
//  of Einstein-Bose quantum gases). For time-discretization one can use either
//  the first-order implicit Euler method or the second-order Crank-Nicolson
//  method.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates.
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi.
//
//  Domain: square (-1, 1)^2.
//
//  BC:  homogeneous Dirichlet everywhere on the boundary.

const int INIT_REF_NUM = 2;                       // Number of initial uniform refinements.
const int P_INIT = 4;                             // Initial polynomial degree.
const double time_step = 0.005;                   // Time step.
const double T_FINAL = 2;                         // Time interval length.
const int TIME_INTEGRATION = 2;                   // 1 for implicit Euler, 2 for Crank-Nicolson.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem constants
const double H = 1;                               // Planck constant 6.626068e-34.
const double M = 1;                               // Mass of boson.
const double G = 1;                               // Coupling constant.
const double OMEGA = 1;                           // Frequency.

// Initial conditions.
scalar init_cond(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-10*(x*x + y*y));
  dx = val * (-20.0 * x);
  dy = val * (-20.0 * y);
  return val;
}

// Boundary markers.
const int BDY_DIRICHLET_1 = 1, BDY_DIRICHLET_2 = 2, BDY_DIRICHLET_3 = 3, BDY_DIRICHLET_4 = 4;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_DIRICHLET_1, BDY_DIRICHLET_2, BDY_DIRICHLET_3, BDY_DIRICHLET_4));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_DIRICHLET_1, BDY_DIRICHLET_2, BDY_DIRICHLET_3, BDY_DIRICHLET_4));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Previous time level solution.
  Solution psi_prev_time(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_INTEGRATION == 1) {
    wf.add_matrix_form(callback(J_euler), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_euler), HERMES_ANY, &psi_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_cranic), HERMES_ANY, &psi_prev_time);
  }

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize views.
  ScalarView view("", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &psi_prev_time, coeff_vec, matrix_solver);
  
  // Time stepping loop:
  int nstep = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    info("Time step %d:", ts);

    // Perform Newton's iteration.
    info("Solving nonlinear problem:");
    bool verbose = true;
    if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
        NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Update previous time level solution.
    Solution::vector_to_solution(coeff_vec, &space, &psi_prev_time);

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time step %d", ts);
    view.set_title(title);
    view.show(&psi_prev_time);
  }

  delete coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
