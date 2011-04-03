#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example shows how to discretize the first-order time-domain Maxwell's equations 
// with vector-valued E (an Hcurl function) and scalar B (an L2 function). Time integration 
// is done using an arbitrary R-K method (see below).
//
// PDE: \partial E / \partial t = SPEED_OF_LIGHT**2 * curl B,
//      \partial B / \partial t = - curl E.
//
// Note: curl E = \partial E_2 / \partial x - \partial E_1 / \partial y
//       curl B = (\partial B / \partial y, - \partial B / \partial x)
//       
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: Rectangular domain (-a, a) x (-b, b)... See mesh file domain.mesh.
//
// BC: perfect conductor for E on the entire boundary.
//
// The following parameters can be changed:

const int P_INIT = 4;                              // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 2;                        // Number of initial uniform mesh refinements.
const double time_step = 0.05;                     // Time step.
const double T_FINAL = 35.0;                       // Final time.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

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
ButcherTableType butcher_table_type = Implicit_RK_1;
//ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;
//ButcherTableType butcher_table_type = Implicit_Radau_IIA_3_5;

// Boundary markers.
const std::string BDY = "Perfect conductor";

// Problem parameters.
const double C_SQUARED = 1;                      // Square of wave speed.                     

// Weak forms.
#include "definitions.cpp"

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
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  Solution B_sln(&mesh, 0.0);
  Hermes::vector<Solution*> slns_time_prev(&E_sln, &B_sln);
  Hermes::vector<Solution*> slns_time_new(&E_sln, &B_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(time_step, C_SQUARED, &E_sln, &B_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace E_space(&mesh, &bcs, P_INIT);
  L2Space B_space(&mesh, P_INIT);

  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&E_space, &B_space)));

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&E_space, &B_space), is_linear);

  // Initialize views.
  ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
  E2_view.fix_scale_width(50);
  ScalarView B_view("Solution B", new WinGeom(0, 405, 400, 350));
  B_view.fix_scale_width(50);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping loop.
  double current_time = time_step; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, slns_time_prev, slns_time_new, false, verbose))
      error("Runge-Kutta time step failed, try to decrease time step size.");

    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time);
    E1_view.set_title(title);
    E1_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time);
    E2_view.set_title(title);
    E2_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);
    sprintf(title, "B, t = %g", current_time);
    B_view.set_title(title);
    B_view.show(&B_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);

    // Update time.
    current_time += time_step;
  
    //View::wait(HERMES_WAIT_KEYPRESS);

  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
