#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

//  This example models a nonstationary distribution of temperature within a wall
//  exposed to ISO fire.
//
//  PDE: non-stationary heat transfer equation
//       HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0.
//  This equation is, however, written in such a way that the time-derivative 
//  is on the left and everything else on the right:
//
//  dT/dt = LAMBDA * Laplace T / (HEATCAP * RHO).
//
//  We only need the weak formulation of the right-hand side.
//
//  Domain: rectangle 4.0 x 0.5 (file wall.mesh).
//
//  IC:  T = TEMP_INIT.
//  BC:  Bottom edge: dT/dn = ALPHA_BOTTOM*(T_fire(x, time) - T)
//       Vertical edges: dT/dn = 0 
//       Top edge: dT/dn = ALPHA_TOP*(TEMP_EXT_TOP - T)
//
//  Time-stepping: Arbitrary Runge-Kutta methods.
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 5;                       // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods: 
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_4_5.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, 
//   Implicit_DIRK_ISMAIL_7_45_embedded. 
ButcherTableType butcher_table_type = Implicit_RK_1;

// Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_RIGHT = 2;
const int BDY_TOP = 3;
const int BDY_LEFT = 4;

// Problem parameters.
const double TEMP_INIT = 20;       // Initial temperature.
const double TEMP_EXT_TOP = 20;    // Exterior temperature top;

const double ALPHA_BOTTOM = 25;    // Heat flux coefficient on the bottom edge.
const double ALPHA_TOP = 8;        // Heat flux coefficient on the top edge.
const double LAMBDA = 1.0;         // Thermal conductivity of the material.
const double HEATCAP = 1020;       // Heat capacity.
const double RHO = 2200;           // Material density.
const double T_FINAL = 4000;       // Length of time interval in seconds.

// Space distribution of fire temperature.
template<typename Real>
Real T_fire_x(Real x) {
  return -1./32 * x*x*x + 3./16 * x*x;
}

// Temporal distribution of fire temperature.
template<typename Real>
Real T_fire_t(Real t) {
  if (0 <= t  &&  t <= 100) return 0;
  if (100 <= t  &&  t <= 600) return 980. / 500 * (t - 100.);
  if (600 <= t  &&  t <= 1800) return 980;
  if (1800 <= t  &&  t <= 3000) return 980 - 980. / 1200 * (t - 1800.);
  return 0.;
}

// Fire temperature as function of x and time.
template<typename Real>
Real T_fire(Real x, Real t) {
  return T_fire_x(x) * T_fire_t(t) + 20;
}

// Weak forms.
#include "forms.cpp"

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
  mloader.load("wall.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_BOTTOM, INIT_REF_NUM_BDY);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_RIGHT, BDY_LEFT));
  bc_types.add_bc_newton(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP));

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, NULL, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Previous and next time level solutions.
  Solution* sln_time_prev = new Solution(&mesh, TEMP_INIT);
  Solution* sln_time_new = new Solution(&mesh);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian_vol));
  wf.add_vector_form(callback(stac_residual_vol), HERMES_ANY, sln_time_prev);
  wf.add_matrix_form_surf(callback(stac_jacobian_bottom), BDY_BOTTOM, sln_time_prev);
  wf.add_vector_form_surf(stac_residual_bottom, stac_residual_bottom_ord, BDY_BOTTOM, sln_time_prev);
  wf.add_matrix_form_surf(callback(stac_jacobian_top), BDY_TOP, sln_time_prev);
  wf.add_vector_form_surf(callback(stac_residual_top), BDY_TOP, sln_time_prev);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Initialize views.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 1000, 250));
  Tview.fix_scale_width(30);

  // Time stepping loop:
  double current_time = time_step; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    bool is_linear = true;
    if (!rk_time_step(current_time, time_step, &bt, sln_time_prev, sln_time_new, &dp, matrix_solver,
		      verbose, is_linear)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    Tview.set_title(title);
    Tview.show(sln_time_new);

    // Copy solution for the new time step.
    sln_time_prev->copy(sln_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  } 
  while (current_time < T_FINAL);

  // Cleanup.
  delete sln_time_prev;
  delete sln_time_new;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
