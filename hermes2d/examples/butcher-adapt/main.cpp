#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"
#include "function/norm.h"

using namespace RefinementSelectors;

//  This example is derived from the example "butcher" and it 
//  shows how adaptive time-stepping can be done with a pair 
//  of embedded Runge-Kutta methods. By embedded we mean that 
//  the Butcher's table has two B rows. After calculating the 
//  stages K_1, K_2, ..., K_s, the two B rows are used to 
//  calculate two different approximations Y_{n+1} on the next
//  time level, with different orders of accuracy. With those
//  one works as usual.  
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
double time_step = 0.2;                            // Time step.
const double T_FINAL = 5.0;                        // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
const double TIME_TOL_LOWER = 0.1;                 // If rel. temporal error is less than this threshold, increase time step
                                                   // but do not repeat time step (this might need further research).
const double TIME_TOL_UPPER = 1.0;                 // If rel. temporal error is greater than this threshold, decrease time 
                                                   // step size and repeat time step.
const double TIME_STEP_INC_RATIO = 1.1;            // Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_DEC_RATIO = 0.5;            // Time step decrease ratio (applied when rel. temporal error is too large).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table:
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, 
// Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, 
// Implicit_Lobatto_IIIC_2_2, Explicit_RK_3, Explicit_RK_4, Implicit_Lobatto_IIIA_3_4, 
// Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_4_5,
// Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, Implicit_DIRK_7_45_embedded. 

ButcherTableType butcher_table_type = Implicit_ESDIRK_TRBDF2_3_23_embedded;
//ButcherTableType butcher_table_type = Implicit_ESDIRK_TRX2_3_23_embedded;
//ButcherTableType butcher_table_type = Implicit_DIRK_7_45_embedded;

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

  // Initialize an error vector for adaptive time stepping.
  scalar* err_vec = new scalar[ndof];

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, space, is_linear);

  // Initialize views.
  OrderView oview("Mesh", new WinGeom(0, 0, 480, 400));
  oview.show(space);
  ScalarView eview("Error", new WinGeom(490, 0, 500, 400));
  ScalarView sview("Solution", new WinGeom(1000, 0, 500, 400));

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    bool is_linear = false;
    if (!rk_time_step(current_time, time_step, &bt, coeff_vec, err_vec, &dp, matrix_solver,
		      verbose, is_linear, NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert err_vec into an error function (Dirichlet lift turned off).
    Solution* error_function = new Solution(&mesh);
    bool add_dir_lift = false;
    Solution::vector_to_solution(err_vec, space, error_function, add_dir_lift);

    // Plot error function.
    // Show the new time level solution.
    char title[100];
    sprintf(title, "Error, t = %g", current_time);
    eview.set_title(title);
    eview.show(error_function, HERMES_EPS_VERYHIGH);

    /* DEBUG - THE B2-ROW IN TABLE "Implicit_DIRK_7_45_embedded" IS PROBABLY
               WRONG, SO ADAPTIVITY IS TEMPORARILY DISABLED.
    // Calculate relative time stepping error and decide whether the 
    // time step can be accepted. If not, then the time step size is 
    // reduced and the entire time step repeated. If yes, then another
    // check is run, and if the relative error is very low, time step 
    // is increased.
    double rel_err = calc_norm(error_function, HERMES_H1_NORM) / calc_norm(sln, HERMES_H1_NORM) * 100;
    if (rel_err > TIME_TOL_UPPER) {
      info("rel_err = %g%% (above upper limit %g%%) -> decreasing time step from %g to %g.", 
           rel_err, TIME_TOL_UPPER, time_step, time_step * TIME_STEP_DEC_RATIO);
      time_step *= TIME_STEP_DEC_RATIO;
      continue;
    }
    if (rel_err < TIME_TOL_LOWER) {
      info("rel_err = %g%% (below lower limit %g%%) -> increasing time step from %g to %g.", 
           rel_err, TIME_TOL_UPPER, time_step, time_step * TIME_STEP_INC_RATIO);
      time_step *= TIME_STEP_INC_RATIO;
      continue;
    }
    */
   
    // Convert coeff_vec into a new time level solution.
    Solution::vector_to_solution(coeff_vec, space, sln);

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
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
