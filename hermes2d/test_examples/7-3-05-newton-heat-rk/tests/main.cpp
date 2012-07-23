#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that example P03-timedep/02-cathedral-rk works correctly.

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 1;                       // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

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

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;         // Thermal conductivity of the material.
const double HEATCAP = 1e2;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 5*time_step;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements();
  mesh.refine_towards_boundary("Boundary_air", INIT_REF_NUM_BDY);
  mesh.refine_towards_boundary("Boundary_ground", INIT_REF_NUM_BDY);

  // Previous and next time level solutions.
  Solution<double>* sln_time_prev = new ConstantSolution<double>(&mesh, TEMP_INIT);
  Solution<double>* sln_time_new = new Solution<double>(&mesh);

  // Initialize the weak formulation.
  double current_time = 0;

  CustomWeakFormHeatRK wf("Boundary_air", ALPHA, LAMBDA, HEATCAP, RHO,
                          &current_time, TEMP_INIT, T_FINAL);

  // Initialize boundary conditions.
 Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("Boundary_ground", TEMP_INIT);
 Hermes::Hermes2D::EssentialBCs<double>bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, &space, &bt);

  runge_kutta.set_freeze_jacobian();
  runge_kutta.set_verbose_output(false);
    
  // Time stepping loop:
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    try
    {
      runge_kutta.setTime(current_time);
      runge_kutta.setTimeStep(time_step);
      runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
    }

    // Copy solution for the new time step.
    sln_time_prev->copy(sln_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
  }
  while (current_time < T_FINAL);

  /* Begin test */

  bool success = true;

  if(fabs(sln_time_new->get_pt_value(-3.5, 17.0) - 10.005262) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(-1.0, 2.0) - 10.0) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(0.0, 9.5) - 9.995515) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value( 1.0, 2.0) - 10.0) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(3.5, 17.0) - 10.005262) > 1E-6) success = false;

  if(success)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }
}