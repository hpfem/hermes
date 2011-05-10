#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example "19-rk-adapt-time-step" works correctly.

const int INIT_GLOB_REF_NUM = 3;                   // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                    // Number of initial refinements towards boundary.
const int P_INIT = 2;                              // Initial polynomial degree.
double time_step = 0.001;                          // Time step.
const double T_FINAL = 3*time_step;                // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
const double TIME_TOL_UPPER = 1.0;                 // If rel. temporal error is greater than this threshold, decrease time 
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
//   Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, Implicit_DIRK_7_45_embedded. 
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "definitions.cpp"

// Initial condition.
#include "initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Check number of command-line parameters.
  if (argc < 2) 
    error("Not enough parameters: Provide a Butcher's table type.");

  int b_type = atoi(argv[1]);
  info ("%d", b_type);

  switch (b_type)
  {
    case 1: butcher_table_type = Explicit_HEUN_EULER_2_12_embedded; break;
    case 2: butcher_table_type = Explicit_BOGACKI_SHAMPINE_4_23_embedded; break;
    case 3: butcher_table_type = Explicit_FEHLBERG_6_45_embedded; break;
    case 4: butcher_table_type = Explicit_CASH_KARP_6_45_embedded; break;
    case 5: butcher_table_type = Explicit_DORMAND_PRINCE_7_45_embedded; break;

    case 6: butcher_table_type = Implicit_ESDIRK_TRBDF2_3_23_embedded; break;
    case 7: butcher_table_type = Implicit_ESDIRK_TRX2_3_23_embedded; break;
    case 8: butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded; break;
    case 9: butcher_table_type = Implicit_SDIRK_CASH_5_24_embedded; break;
    case 10: butcher_table_type = Implicit_SDIRK_CASH_5_34_embedded; break;
    case 11: butcher_table_type = Implicit_DIRK_ISMAIL_7_45_embedded; break;

    default: error("Admissible command-line options are from 1 to 11.");
  }

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

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
    bool verbose = false;
    bool jacobian_changed = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, &sln_time_prev, 
                                  &sln_time_new, &time_error_fn, jacobian_changed, verbose, 
                                  NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

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

    // Copy solution for next time step.
    sln_time_prev.copy(&sln_time_new);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  info("Coordinate (-8.0, -8.0) value = %lf", sln_time_prev.get_pt_value(-8.0, -8.0));
  info("Coordinate (-5.0, -5.0) value = %lf", sln_time_prev.get_pt_value(-5.0, -5.0));
  info("Coordinate (-3.0, -3.0) value = %lf", sln_time_prev.get_pt_value(-3.0, -3.0));
  info("Coordinate ( 0.0,  0.0) value = %lf", sln_time_prev.get_pt_value(0.0,  0.0));
  info("Coordinate ( 3.0,  3.0) value = %lf", sln_time_prev.get_pt_value(3.0,  3.0));
  info("Coordinate ( 5.0,  5.0) value = %lf", sln_time_prev.get_pt_value(5.0,  5.0));
  info("Coordinate ( 8.0,  8.0) value = %lf", sln_time_prev.get_pt_value(8.0,  8.0));

  double coor_x_y[7] = {-8.0, -5.0, -3.0, 0.0, 3.0, 5.0, 8.0};
  bool success = true;

  switch (b_type)
  {
    case 1: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695536) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.260018) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279160) > 1E-6) success = false;
            break;

    case 2: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003562) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259981) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279212) > 1E-6) success = false;
            break;

    case 3: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259990) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279225) > 1E-6) success = false;
            break;

    case 4: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259989) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279225) > 1E-6) success = false;
            break;

    case 5: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259989) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279226) > 1E-6) success = false;
            break;

    case 6: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003562) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695539) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259981) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279250) > 1E-6) success = false;
            break;

    case 7: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003562) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695539) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259986) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279238) > 1E-6) success = false;
            break;

    case 8: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259987) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279225) > 1E-6) success = false;
            break;

    case 9: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259989) > 1E-6) success = false;
            if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279224) > 1E-6) success = false;
            break;

    case 10: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259989) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279224) > 1E-6) success = false;
             break;

    case 11: if (fabs(sln_time_prev.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.044235) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.253124) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.493350) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.003563) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.695538) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.259989) > 1E-6) success = false;
             if (fabs(sln_time_prev.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.279226) > 1E-6) success = false;
             break;

    default: error("Admissible command-line options are from 1 to 11.");
  }

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
