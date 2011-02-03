#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

using namespace RefinementSelectors;

//  This example is derived from the tutorial example 19.
//  It uses general Butcher's tables to perform arbitrary 
//  explicit or implicit low-order or higher-order Runge-Kutta
//  time integration. Example 19 can just do implicit Euler.
//  The model problem is a simple nonlinear parabolic PDE.
//
//  The function rk_time_step() needs more optimisation, see
//  a todo list at the beginning of file src/runge-kutta.cpp.
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

const int INIT_GLOB_REF_NUM = 3;                   // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;                    // Number of initial refinements towards boundary.
const int P_INIT = 2;                              // Initial polynomial degree.
const double time_step = 0.2;                      // Time step.
const double T_FINAL = 5.0;                        // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, 
// Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, 
// Implicit_Lobatto_IIIC_2_2, Explicit_RK_3, Explicit_RK_4, Implicit_Lobatto_IIIA_3_4, 
// Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_4_5,
// Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
// Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, Implicit_DIRK_7_45_embedded. 

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u) 
{ 
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) 
{ 
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
{ return dir_lift(x, y, dx, dy);}

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
Real heat_src(Real x, Real y) { return 1.0;}

// Weak forms.
#include "forms.cpp"

// Main function.
int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

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
  Solution* sln = new Solution(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian));
  wf.add_vector_form(callback(stac_residual));

  // Project the initial condition on the FE space to obtain initial solution coefficient vector.
  info("Projecting initial condition to translate initial condition into a vector.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(space, sln, coeff_vec, matrix_solver);

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
    bool is_linear = false;
    if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver,
		      verbose, is_linear, NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert coeff_vec into a new time level solution.
    Solution::vector_to_solution(coeff_vec, space, sln);

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(sln, HERMES_EPS_VERYHIGH);
    oview.show(space);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  // Cleanup.
  delete [] coeff_vec;
  delete space;
  delete sln;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
