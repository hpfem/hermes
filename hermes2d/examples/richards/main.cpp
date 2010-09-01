#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

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

const int INIT_GLOB_REF_NUM = 3;       // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;        // Number of initial refinements towards boundary.
const int P_INIT = 3;                  // Initial polynomial degree.
const double TAU = 5e-3;               // Time step.
const double T_FINAL = 0.4;            // Time interval length.
const int TIME_INTEGRATION = 1;        // 1... implicit Euler, 2... Crank-Nicolson.
const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;       // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

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

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = (100 - 2*x)/2.5 * pow(y/100, Y_POWER);
  dy = x*(100 - x)/2.5 * pow(y/100, Y_POWER - 1) * 1./100;
  return x*(100 - x)/2.5 * pow(y/100, Y_POWER) - 1000;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return init_cond(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Solution for the previous time level.
  Solution u_prev_time;

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, H2D_UNSYM, H2D_ANY, &u_prev_time);
    wf.add_vector_form(res_euler, res_ord, H2D_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, H2D_UNSYM, H2D_ANY, &u_prev_time);
    wf.add_vector_form(res_cranic, res_ord, H2D_ANY, &u_prev_time);
  }

  // Initialize matrix solver.
  Matrix* mat; Vector* coeff_vec; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, coeff_vec, solver);

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  Solution* sln_tmp = new Solution(&mesh, init_cond);
  project_global(&space, H2D_H1_NORM, sln_tmp, &u_prev_time, coeff_vec);
  delete sln_tmp;

  // Initialize views.
  ScalarView sview("Solution", 0, 0, 500, 400);
  OrderView oview("Mesh", 520, 0, 450, 400);
  oview.show(&space);
  sview.show(&u_prev_time);
  //View::wait(H2DV_WAIT_KEYPRESS);

  // Time stepping loop:
  double current_time = 0.0;
  int t_step = 1;
  do {
    info("---- Time step %d, t = %g s.", t_step, current_time); t_step++;

    // Newton's method.
    info("Performing Newton's method.");
    bool verbose = true; // Default is false.
    if (!solve_newton(&space, &wf, coeff_vec, matrix_solver, NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
      error("Newton's method did not converge.");

    // Update previous time level solution.
    u_prev_time.set_fe_solution(&space, coeff_vec);

    // Update time.
    current_time += TAU;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&u_prev_time);
  } while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

