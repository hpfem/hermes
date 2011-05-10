#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example is analogous to example 01-timedep-adapt-space-only except that 
//  a fixed mesh is used and only the time stepping is adaptive. An arbitrary 
//  embedded Runge-Kutta method can be used. By embedded we mean that the 
//  Butcher's table contains two B rows. After calculating the stages K_1, 
//  K_2, ..., K_s, the two B rows are used to calculate two different 
//  approximations Y_{n+1} on the next time level, with different orders 
//  of accuracy. With those one works as usual.
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
double time_step = 0.05;                           // Time step.
const double T_FINAL = 5.0;                        // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
const double TIME_TOL_UPPER = 5.0;                 // If rel. temporal error is greater than this threshold, decrease time 
                                                   // step size and repeat time step.
const double TIME_TOL_LOWER = 0.5;                 // If rel. temporal error is less than this threshold, increase time step
                                                   // but do not repeat time step (this might need further research).
const double TIME_STEP_INC_RATIO = 1.1;            // Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_DEC_RATIO = 0.8;            // Time step decrease ratio (applied when rel. temporal error is too large).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double ALPHA = 4.0;                          // For the nonlinear thermal conductivity.

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, 
//   Implicit_DIRK_ISMAIL_7_45_embedded. 
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

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

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential(BDY_DIRICHLET);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Convert initial condition into a Solution.
  InitialSolutionHeatTransfer sln_time_prev(&mesh);
  Solution sln_time_new(&mesh);
  Solution time_error_fn(&mesh, 0.0);

  // Initialize the weak formulation
  WeakFormHeatTransferNewtonTimedep wf(ALPHA, time_step, &sln_time_prev);

  // Initialize the discrete problem.
  DiscreteProblem dp(&wf, &space);

  // Initialize views.
  ScalarView sview_high("Solution (higher-order)", new WinGeom(0, 0, 500, 400));
  ScalarView sview_low("Solution (lower-order)", new WinGeom(490, 0, 500, 400));
  ScalarView eview("Temporal error", new WinGeom(1000, 0, 500, 400));
  eview.fix_scale_width(50);

  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  info("Time step history will be saved to file time_step_history.dat.");

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    bool jacobian_changed = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, &sln_time_prev, 
                                  &sln_time_new, &time_error_fn, jacobian_changed, verbose, 
                                  NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Plot error function.
    char title[100];
    sprintf(title, "Temporal error, t = %g", current_time + time_step);
    eview.set_title(title);
    AbsFilter abs_tef(&time_error_fn);
    eview.show(&abs_tef, HERMES_EPS_VERYHIGH);

    // Calculate relative time stepping error and decide whether the 
    // time step can be accepted. If not, then the time step size is 
    // reduced and the entire time step repeated. If yes, then another
    // check is run, and if the relative error is very low, time step 
    // is increased.
    double rel_err_time = hermes2d.calc_norm(&time_error_fn, HERMES_H1_NORM) / 
                          hermes2d.calc_norm(&sln_time_new, HERMES_H1_NORM) * 100;
    info("rel_err_time = %g%%", rel_err_time);
    if (rel_err_time > TIME_TOL_UPPER) {
      info("rel_err_time above upper limit %g%% -> decreasing time step from %g to %g and repeating time step.", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_DEC_RATIO);
      time_step *= TIME_STEP_DEC_RATIO;
      continue;
    }
    if (rel_err_time < TIME_TOL_LOWER) {
      info("rel_err_time = below lower limit %g%% -> increasing time step from %g to %g", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_INC_RATIO);
      time_step *= TIME_STEP_INC_RATIO;
    }
   
    // Add entry to the timestep graph.
    time_step_graph.add_values(current_time, time_step);
    time_step_graph.save("time_step_history.dat");

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    sprintf(title, "Solution (higher-order), t = %g", current_time);
    sview_high.set_title(title);
    sview_high.show(&sln_time_new, HERMES_EPS_HIGH);
    sprintf(title, "Solution (lower-order), t = %g", current_time);
    sview_low.set_title(title);
    SumFilter sln_time_new_low(Hermes::vector<MeshFunction*>(&sln_time_new, &time_error_fn), 
                               Hermes::vector<int>(H2D_FN_VAL, H2D_FN_VAL));
    sview_low.show(&sln_time_new_low, HERMES_EPS_HIGH);

    // Copy solution for next time step.
    sln_time_prev.copy(&sln_time_new);

    // Increase counter of time steps.
    ts++;

    //View::wait(HERMES_WAIT_KEYPRESS);

  } 
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
