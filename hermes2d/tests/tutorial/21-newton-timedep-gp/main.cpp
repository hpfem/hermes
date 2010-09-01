#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
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
//  BC:  homogeneous Dirichlet everywhere on the boundary

const int INIT_REF_NUM = 2;      // Number of initial uniform refinements.
const int P_INIT = 4;            // Initial polynomial degree.
const double TAU = 0.005;        // Time step.
const double T_FINAL = 2*TAU + 1e-4;        // Time interval length.
const int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson.
const double NEWTON_TOL = 1e-5;  // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100; // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem constants
const double H = 1;              // Planck constant 6.626068e-34.
const double M = 1;              // Mass of boson.
const double G = 1;              // Coupling constant.
const double OMEGA = 1;          // Frequency.

// Initial conditions.
scalar init_cond(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-10*(x*x + y*y));
  dx = val * (-20.0 * x);
  dy = val * (-20.0 * y);
  return val;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
 return 0;
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
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Previous time level solution.
  Solution psi_prev_time(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(J_euler), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_euler), H2D_ANY, &psi_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_cranic), H2D_ANY, &psi_prev_time);
  }

  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  bool is_complex = true;
  Vector* coeff_vec = new AVector(0, is_complex); 
  project_global(space, H2D_H1_NORM, &psi_prev_time, Tuple<Solution*>(), coeff_vec, is_complex);

  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {

    info("Time step %d:", ts);

    // Newton's method.
    info("Performing Newton's method.");
    bool verbose = true; // Default is false.
    if (!solve_newton(space, &wf, coeff_vec, matrix_solver, 
		      NEWTON_TOL, NEWTON_MAX_ITER, verbose, is_complex))
      error("Newton's method did not converge.");

    // Update previous time level solution.
    psi_prev_time.set_fe_solution(space, coeff_vec);
  }

  delete coeff_vec;

  #define ERROR_SUCCESS                                0
  #define ERROR_FAILURE                               -1
  printf("Success!\n");
  return ERROR_SUCCESS;
}
