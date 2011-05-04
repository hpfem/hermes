#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that the example "maxwell/resonator-time-domain-I" works correctly.

// The following parameters can be changed:
const int P_INIT = 6;                              // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 1;                        // Number of initial uniform mesh refinements.
const double time_step = 0.05;                     // Time step.
const double T_FINAL = time_step*5;                // Final time.
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
#include "../definitions.cpp"

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
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  Solution B_sln(&mesh, 0.0);
  Hermes::vector<Solution*> slns(&E_sln, &B_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(C_SQUARED);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY, 0.0);
  EssentialBCs bcs_E(&bc_essential);
  EssentialBCs bcs_B;

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace E_space(&mesh, &bcs_E, P_INIT);
  H1Space B_space(&mesh, &bcs_B, P_INIT);
  //L2Space B_space(&mesh, P_INIT);
  Hermes::vector<Space *> spaces = Hermes::vector<Space *>(&E_space, &B_space);
  info("ndof = %d.", Space::get_num_dofs(spaces));

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, spaces);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping loop.
  double current_time = time_step; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, slns, slns, jacobian_changed, verbose))
      error("Runge-Kutta time step failed, try to decrease time step size.");

    // Update time.
    current_time += time_step;
  
  } while (current_time < T_FINAL);

  double coord_x[4] = {0.3, 0.6, 0.9, 1.4};
  double coord_y[4] = {0, 0.3, 0.5, 0.7};

  info("Coordinate (0.3, 0.0) value = %lf", B_sln.get_pt_value(coord_x[0], coord_y[0]));
  info("Coordinate (0.6, 0.3) value = %lf", B_sln.get_pt_value(coord_x[1], coord_y[1]));
  info("Coordinate (0.9, 0.5) value = %lf", B_sln.get_pt_value(coord_x[2], coord_y[2]));
  info("Coordinate (1.4, 0.7) value = %lf", B_sln.get_pt_value(coord_x[3], coord_y[3]));

  double t_value[4] = {0.0, -0.065100, -0.146515, -0.247677};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (fabs(t_value[i] - B_sln.get_pt_value(coord_x[i], coord_y[i])) > 1E-6) success = false;
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
