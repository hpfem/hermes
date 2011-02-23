#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 22-newton-timedep-gp works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform refinements.
const int P_INIT = 4;                             // Initial polynomial degree.
const double TAU = 0.005;                         // Time step.
const double T_FINAL = 2*TAU + 1e-4;              // Time interval length.
const int TIME_DISCR = 2;                         // 1 for implicit Euler, 2 for Crank-Nicolson.
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
  if(TIME_DISCR == 1) {
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

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &psi_prev_time, coeff_vec, matrix_solver);
  
  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
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
  }

  delete coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  AbsFilter mag2(&psi_prev_time);
  int success = 1;
  double eps = 1e-5;
  double val = std::abs(mag2.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.1, 0.1) psi value = %lf", val);
  if (fabs(val - (0.804900)) > eps) {
    printf("Coordinate ( 0.1, 0.1) psi value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.1, -0.1));
  info("Coordinate ( 0.1, -0.1) psi value = %lf", val);
  if (fabs(val - (0.804900)) > eps) {
    printf("Coordinate ( 0.1, -0.1) psi value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.2, 0.1));
  info("Coordinate ( 0.2, 0.1) psi value = %lf", val);
  if (fabs(val - (0.602930)) > eps) {
    printf("Coordinate ( 0.2, 0.1) psi value = %lf\n", val);
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
