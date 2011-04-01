#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves a time-domain resonator problem for the Maxwell's equation. 
// It is very similar to resonator-time-domain-I but B is eliminated from the 
// equations, thus converting the first-order system into one second -order
// equation in time. The second-order equation in time is decomposed back into 
// a first-order system in time in the standard way (see example wave-1). Time 
// discretization is performed using the implicit Euler method.
//
// The function rk_time_step() needs more optimisation, see a todo list at 
// the beginning of file src/runge-kutta.h.
//
// PDE: \frac{1}{SPEED_OF_LIGHT**2}\frac{\partial^2 E}{\partial t^2} + curl curl E = 0,
// converted into
//
//      \frac{\partial E}{\partial t} = F,
//      \frac{\partial F}{\partial t} = - SPEED_OF_LIGHT**2 * curl curl E.
//
// Approximated by
// 
//      \frac{E^{n+1}}{tau}                   - F^{n+1}             = \frac{E^{n}}{tau},
//      SPEED_OF_LIGHT**2 * curl curl E^{n+1} + \frac{F^{n+1}}{tau} = \frac{F^{n}}{tau}.
//
// Domain: Rectangular domain (-pi/2, pi/2) x (-pi/2, pi/2)... See mesh file domain.mesh.
//
// BC:  E \times \nu = 0 on the boundary (perfect conductor),
//      F \times \nu = 0 on the boundary (E \times \nu = 0 => \partial E / \partial t \times \nu = 0).
//
// IC:  Prescribed wave for E, zero for F.
//
// The following parameters can be changed:

const int P_INIT = 2;                              // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 4;                        // Number of initial uniform mesh refinements.
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
const double C_SQUARED = 1;                      // Square of wave speed.                     

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  Solution F_sln(&mesh, 0.0, 0.0);
  Hermes::vector<Solution*> slns_time_prev(&E_sln, &F_sln);
  Hermes::vector<Solution*> slns_time_new(&E_sln, &F_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(time_step, C_SQUARED, &E_sln, &F_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace E_space(&mesh, &bcs, P_INIT);
  HcurlSpace F_space(&mesh, &bcs, P_INIT);

  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&E_space, &F_space)));

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&E_space, &F_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Initialize views.
  ScalarView E1_view("Solution E1", new WinGeom(0, 0, 420, 300));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("Solution E2", new WinGeom(430, 0, 420, 300));
  E2_view.fix_scale_width(50);
  ScalarView F1_view("Solution F1", new WinGeom(0, 355, 420, 300));
  F1_view.fix_scale_width(50);
  ScalarView F2_view("Solution E2", new WinGeom(430, 355, 420, 300));
  F2_view.fix_scale_width(50);

  // Time stepping loop.
  double current_time = 0; int ts = 1; bool rhs_only = false;
  do
  {
    // Perform one implicit Euler time step.
    info("Implicit Euler time step (t = %g s, time_step = %g s).", current_time, time_step);
    bool verbose = true;

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
    rhs_only = true;

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space*>(&E_space, &F_space), 
                                                      Hermes::vector<Solution*>(&E_sln, &F_sln));
    else error ("Matrix solver failed.\n");

    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time);
    E1_view.set_title(title);
    E1_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time);
    E2_view.set_title(title);
    E2_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);
    sprintf(title, "F1, t = %g", current_time);
    F1_view.set_title(title);
    F1_view.show(&F_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "F2, t = %g", current_time);
    F2_view.set_title(title);
    F2_view.show(&F_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);

    // Update time.
    current_time += time_step;
  
    //View::wait(HERMES_WAIT_KEYPRESS);

  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
