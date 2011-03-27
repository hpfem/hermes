#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

//  This example solves a simple version of the time-dependent
//  Richard's equation using the backward Euler method in time 
//  combined with the Newton's method in each time step.
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: square (0, 100)^2.
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
//#define CONSTITUTIVE_GENUCHTEN

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 3;                             // Initial polynomial degree.
const double TAU = 5e-3;                          // Time step.
const double T_FINAL = 0.4;                       // Time interval length.
const int TIME_INTEGRATION = 1;                   // 1... implicit Euler, 2... Crank-Nicolson.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// For the definition of initial condition.
int Y_POWER = 10;

// Problem parameters.
double K_S = 20.464;
double ALPHA = 0.001;
double THETA_R = 0;
double THETA_S = 0.45;
double H_R = -1000;
double A = 100;
double L = 100;
double STORATIVITY = 0.0;
double M = 0.6;
double N = 2.5;

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = (100 - 2*x)/2.5 * pow(y/100, Y_POWER);
  dy = x*(100 - x)/2.5 * pow(y/100, Y_POWER - 1) * 1./100;
  return x*(100 - x)/2.5 * pow(y/100, Y_POWER) - 1000;
}

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(double x, double y)
{
  double dx, dy;
  return init_cond(x, y, dx, dy);
}

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Solution for the previous time level.
  Solution u_prev_time;

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, HERMES_NONSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_euler, res_ord, HERMES_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, HERMES_NONSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_cranic, res_ord, HERMES_ANY, &u_prev_time);
  }

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)];
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  OGProjection::project_global(&space, init_cond, coeff_vec, matrix_solver);
  Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(520, 0, 450, 400));
  oview.show(&space);
  sview.show(&u_prev_time);
  
  // Time stepping loop:
  double current_time = 0.0;
  int t_step = 1;
  do {
    info("---- Time step %d, t = %g s.", t_step, current_time); t_step++;

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, &space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Perform Newton's iteration.
    bool verbose = true;
    if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
        NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

    // Cleanup.
    delete matrix;
    delete rhs;
    delete solver;

    // Update time.
    current_time += TAU;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&u_prev_time);
  } while (current_time < T_FINAL);
  
  delete [] coeff_vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

