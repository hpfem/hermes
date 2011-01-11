#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method. The St. Vitus Cathedral
//  in Prague (http://en.wikipedia.org/wiki/St._Vitus_Cathedral)
//  responds to changes in the surrounding air temperature
//  during one 24-hour cycle. You will also learn how to use the
//  solution calculated in the previous time step.
//
//  PDE: non-stationary heat transfer equation
//  HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0.
//
//  Domain: St. Vitus cathedral (cathedral.mesh).
//
//  IC:  T = TEMP_INIT.
//  BC:  T = TEMP_INIT on the bottom edge ... Dirichlet,
//       dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent.
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

const int P_INIT = 4;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 3e2;                     // Time step in seconds.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table:
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2, Implicit_SDIRK_2, 
// Implicit_Lobatto_IIIA_2, Implicit_Lobatto_IIIB_2, Implicit_Lobatto_IIIC_2, Explicit_RK_3, Explicit_RK_4,
// Implicit_Lobatto_IIIA_4, Implicit_Lobatto_IIIB_4, Implicit_Lobatto_IIIC_4. 
ButcherTableType butcher_table_type = Implicit_RK_1;

// Boundary markers.
const std::string BDY_GROUND = "Boundary ground";
const std::string BDY_AIR = "Boundary air";

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e5;         // Thermal conductivity of the material.
const double HEATCAP = 1e6;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 86400;      // Length of time interval (24 hours) in seconds.

// Time-dependent exterior temperature.
template<typename Real>
Real temp_ext(Real t) {
  return TEMP_INIT + 10. * sin(2*M_PI*t/T_FINAL);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Choose a Butcher's table.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_AIR, INIT_REF_NUM_BDY);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_GROUND);
  bc_types.add_bc_newton(BDY_AIR);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_const(BDY_GROUND, TEMP_BND);

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Initialize the solution.
  Solution u_prev_time;


  // Set the initial condition.
  u_prev_time.set_const(&mesh, TEMP_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, BDY_AIR);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, BDY_AIR);

  // Project the initial condition on the FE space to obtain initial solution coefficient vector.
  info("Projecting initial condition to translate initial condition into a vector.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Initialize views.
  double current_time = 0.0; 
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", current_time, temp_ext(current_time));
  //Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // Time stepping loop:
  int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver,
		      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert coeff_vec into a new time level solution.
    Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

    // Update time.
    current_time += time_step;

    // Visualize the solution.
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", current_time, temp_ext(current_time));
    Tview.set_title(title);
    Tview.show(&u_prev_time);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  // Cleanup.
  delete [] coeff_vec;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
