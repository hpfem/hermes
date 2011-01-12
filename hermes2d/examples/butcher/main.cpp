#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

using namespace RefinementSelectors;

//  This example uses general Butcher's tables to perform 
//  arbitrary explicit or implicit low-order or higher-order
//  time integration. The model problem is a simple nonlinear 
//  parabolic PDE.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f.
//
//  Domain: square (-10,10)^2.
//
//  BC: Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int INIT_GLOB_REF_NUM = 4;                   // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 3;                    // Number of initial refinements towards boundary.
const int P_INIT = 2;                              // Initial polynomial degree.
const double time_step = 0.2;                      // Time step.
const double T_FINAL = 5.0;                        // Time interval length.
const double NEWTON_TOL = 1e-3;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table:
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2, Implicit_SDIRK_2, 
// Implicit_Lobatto_IIIA_2, Implicit_Lobatto_IIIB_2, Implicit_Lobatto_IIIC_2, Explicit_RK_3, Explicit_RK_4,
// Implicit_Lobatto_IIIA_4, Implicit_Lobatto_IIIB_4, Implicit_Lobatto_IIIC_4. 

//ButcherTableType butcher_table_type = Implicit_SDIRK_2;
//ButcherTableType butcher_table_type = Implicit_Crank_Nicolson_2;
ButcherTableType butcher_table_type = Implicit_RK_1;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

// Main function.
int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);   

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  Solution u_prev_time(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian), HERMES_NONSYM, HERMES_ANY);
  wf.add_vector_form(callback(stac_residual), HERMES_ANY);

  // Project the initial condition on the FE space to obtain initial solution coefficient vector.
  info("Projecting initial condition to translate initial condition into a vector.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(space, &u_prev_time, coeff_vec, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, space, is_linear);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(510, 0, 460, 400));
  oview.show(space);

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver,
		      verbose, NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert coeff_vec into a new time level solution.
    Solution::vector_to_solution(coeff_vec, space, &u_prev_time);

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&u_prev_time);
    oview.show(space);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  // Cleanup.
  delete [] coeff_vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
