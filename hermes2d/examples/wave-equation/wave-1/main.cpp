#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves a simple linear wave equation by converting it 
// into a system of two first-order equations in time. Time discretization 
// is performed using arbitrary (explicit or implicit, low-order or higher-order) 
// Runge-Kutta methods entered via their Butcher's tables.
// For a list of available R-K methods see the file hermes_common/tables.h.
//
// The function rk_time_step() needs more optimisation, see a todo list at 
// the beginning of file src/runge-kutta.h.
//
// PDE: \frac{1}{C_SQUARED}\frac{\partial^2 u}{\partial t^2} - \Delta u = 0,
// converted into
//
//      \frac{\partial u}{\partial t} = v,
//      \frac{\partial v}{\partial t} = C_SQUARED * \Delta u.
//
// BC:  u = 0 on the boundary,
//      v = 0 on the boundary (u = 0 => \partial u / \partial t = 0).
//
// IC:  smooth peak for u, zero for v.
//
// The following parameters can be changed:

const int P_INIT = 6;                              // Initial polynomial degree of all elements.
const double time_step = 0.01;                     // Time step.
const double T_FINAL = 2.15;                       // Final time.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,

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

// Boundary markers.
const std::string BDY = "1";

// Problem parameters.
const double C_SQUARED = 100;                      // Square of wave speed.                     

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

  // Convert to quadrilaterals.
  //mesh.convert_triangles_to_quads();

  // Refine towards boundary.
  mesh.refine_towards_boundary(BDY, 1, true, true);

  // Refine once towards vertex #4.
  mesh.refine_towards_vertex(4, 1);

  // Initialize solutions.
  CustomInitialConditionWave u_sln(&mesh);
  Solution v_sln(&mesh, 0.0);
  Hermes::vector<Solution*> slns_time_prev(&u_sln, &v_sln);
  Hermes::vector<Solution*> slns_time_new(&u_sln, &v_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(time_step, C_SQUARED, &u_sln, &v_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u_space(&mesh, &bcs, P_INIT);
  H1Space v_space(&mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)));

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u_space, &v_space), is_linear);

  // Initialize views.
  ScalarView u_view("Solution u", new WinGeom(0, 0, 500, 400));
  //u_view.show_mesh(false);
  u_view.fix_scale_width(50);
  ScalarView v_view("Solution v", new WinGeom(510, 0, 500, 400));
  //v_view.show_mesh(false);
  v_view.fix_scale_width(50);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);


  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;

    if (!runge_kutta.rk_time_step(current_time, time_step, slns_time_prev, 
                                  slns_time_new, jacobian_changed, verbose))
      error("Runge-Kutta time step failed, try to decrease time step size.");

    // Visualize the solutions.
    char title[100];
    sprintf(title, "Solution u, t = %g", current_time);
    u_view.set_title(title);
    u_view.show(&u_sln);
    sprintf(title, "Solution v, t = %g", current_time);
    v_view.set_title(title);
    v_view.show(&v_sln);

    // Update time.
    current_time += time_step;

  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
