#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example solves the time-dependent Richard's equation using 
//  adaptive time integration (no dynamical meshes in space yet).
//  Many different time stepping methods can be used. The nonlinear 
//  solver in each time step is the Newton's method. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: rectangle (0, 8) x (0, 6.5).
//  Units: length: cm
//         time: days
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// Constitutive relations.
#define CONSTITUTIVE_GENUCHTEN                    // Van Genuchten or Gardner.

// Choose full domain or half domain.
//const char* mesh_file = "domain-full.mesh";
const char* mesh_file = "domain-half.mesh";

// Adaptive time stepping.
double time_step = 0.3;                           // Time step (in days).
const double time_tol_upper = 1.0;                // If rel. temporal error is greater than this threshold, decrease time 
                                                  // step size and repeat time step.
const double time_tol_lower = 0.5;                // If rel. temporal error is less than this threshold, increase time step
                                                  // but do not repeat time step (this might need further research).
double time_step_dec = 0.8;                       // Timestep decrease ratio after unsuccessful nonlinear solve.
double time_step_inc = 1.1;                       // Timestep increase ratio after successful nonlinear solve.
double time_step_min = 1e-8; 			  // Computation will stop if time step drops below this value. 
                       
// Elements orders and initial refinements.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY_TOP = 1;               // Number of initial mesh refinements towards the top edge.

// Matrix solver.
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
//   Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, Implicit_DIRK_7_45_embedded. 
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

// Newton's method.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for Newton on fine mesh.
int NEWTON_MAX_ITER = 10;                         // Maximum allowed number of Newton iterations.

// Times.
const double STARTUP_TIME = 5.0;                  // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 1000.0;                    // Time interval length.
const double PULSE_END_TIME = 1000.0;             // Time interval of the top layer infiltration.
double current_time = time_step;                  // Global time variable initialized with first time step.

// Problem parameters.
double H_INIT = -15.0;                            // Initial pressure head.
double H_ELEVATION = 10.0;                        // Top constant pressure head -- an infiltration experiment.
const double K_S_vals[4] = {350.2, 712.8, 1.68, 18.64}; 
const double ALPHA_vals[4] = {0.01, 1.0, 0.01, 0.01};
const double N_vals[4] = {2.5, 2.0, 1.23, 2.5};
const double M_vals[4] = {0.864, 0.626, 0.187, 0.864};

const double THETA_R_vals[4] = {0.064, 0.0, 0.089, 0.064};
const double THETA_S_vals[4] = {0.14, 0.43, 0.43, 0.24};
const double STORATIVITY_vals[4] = {0.1, 0.1, 0.1, 0.1};


// Precalculation of constitutive tables.
const int MATERIAL_COUNT = 4;
const int CONSTITUTIVE_TABLE_METHOD = 2;          // 0 - constitutive functions are evaluated directly (slow performance).
					          // 1 - constitutive functions are linearly approximated on interval 
                                                  //     <TABLE_LIMIT; LOW_LIMIT> (very efficient CPU utilization less 
                                                  //     efficient memory consumption (depending on TABLE_PRECISION)).
						  // 2 - constitutive functions are aproximated by quintic splines.
						  
//!Use only if 	CONSTITUTIVE_TABLE_METHOD == 2 !//					  
const int NUM_OF_INTERVALS = 16;                  // Number of intervals.                      
const double INTERVALS_4_APPROX[16] = 
      {-1.0, -2.0, -3.0, -4.0, -5.0, -8.0, -10.0, -12.0, // Low limits of intervals approximated by quintic splines.
      -15.0, -20.0, -30.0, -50.0, -75.0, -100.0,-300.0, -1000.0}; 
int* POL_SEARCH_HELP;                             // This array contains for each integer of h function appropriate 
                                                  // polynomial ID.
double**** K_POLS;                                // First DIM is the interval ID, second DIM is the material ID, 
                                                  // third DIM is the derivative degree, fourth DIM are the coefficients.
double**** C_POLS;
//!------------------------------------------!//


//!Use only if CONSTITUTIVE_TABLE_METHOD == 1 !//
double TABLE_LIMIT = -1000.0; 		          // Limit of precalculated functions (should be always negative value lower 
						  // then the lowest expect value of the solution (consider DMP!!)
const double TABLE_PRECISION = 0.1;               // Precision of precalculated table use 1.0, 0,1, 0.01, etc.....
double** K_TABLE;                                  
double** dKdh_TABLE;
double** ddKdhh_TABLE;
double** C_TABLE;
double** dCdh_TABLE;
bool CONSTITUTIVE_TABLES_READY = false;
double*** POLYNOMIALS;                            // Polynomial approximation of the K(h) function close to saturation.
                                                  // This function has singularity in its second derivative.
						  // First dimension is material ID
						  // Second dimension is the polynomial derivative.
						  // Third dimension are the polynomial's coefficients.
						  
						  
const double LOW_LIMIT = -1.0;                    // Lower bound of K(h) function approximated by polynomials.
const int NUM_OF_INSIDE_PTS = 0;
//!------------------------------------------!//


bool POLYNOMIALS_READY = false;
bool POLYNOMIALS_ALLOCATED = false;


// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M, STORATIVITY;

// Choose here which constitutive relations should be used.
#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
const int BDY_TOP = 1;
const int BDY_RIGHT = 2;
const int BDY_BOTTOM = 3;
const int BDY_LEFT = 4;

// Initial condition.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = 0;
  dy = 0;
  return H_INIT;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y, double time)
{
  if (time < STARTUP_TIME)
    return H_INIT + time/STARTUP_TIME*(H_ELEVATION-H_INIT);
  else if (time > PULSE_END_TIME)
    return H_INIT;
  else
    return H_ELEVATION;
}

// Weak forms.
#include "forms.cpp"

// Additional functionality.
#include "extras.cpp"

// Main function.
int main(int argc, char* argv[])
{
  // Either use exact constitutive relations (slow) (method 0) or precalculate 
  // their linear approximations (faster) (method 1) or
  // precalculate their quintic polynomial approximations (method 2) -- managed by 
  // the following loop "Initializing polynomial approximation".
  if (CONSTITUTIVE_TABLE_METHOD == 1)
    CONSTITUTIVE_TABLES_READY = get_constitutive_tables(1);  // 1 stands for the Newton's method.
  // Points to be used for polynomial approximation of K(h).
  double* points = new double[NUM_OF_INSIDE_PTS];

  // The van Genuchten + Mualem K(h) function is approximated by polynomials close 
  // to zero in case of CONSTITUTIVE_TABLE_METHOD==1.
  // In case of CONSTITUTIVE_TABLE_METHOD==2, all constitutive functions are approximated by polynomials.
  info("Initializing polynomial approximations.");
  for (int i=0; i < MATERIAL_COUNT; i++) {
    init_polynomials(6 + NUM_OF_INSIDE_PTS, LOW_LIMIT, points, NUM_OF_INSIDE_PTS, i);
  }
  POLYNOMIALS_READY = true;
  if (CONSTITUTIVE_TABLE_METHOD == 2) {
    CONSTITUTIVE_TABLES_READY = true ;
    //Assign table limit to global definition.
    TABLE_LIMIT = INTERVALS_4_APPROX[NUM_OF_INTERVALS-1];
  }
  
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load(mesh_file, &basemesh);
  
  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_TOP, INIT_REF_NUM_BDY_TOP);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_TOP);
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_RIGHT, BDY_BOTTOM, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values(&current_time);
  bc_values.add_timedep_function(BDY_TOP, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Solution (initialized by the initial condition) and error function.
  Solution* sln_time_prev = new Solution(&mesh, init_cond);
  Solution* sln_time_new = new Solution(&mesh);
  Solution* time_error_fn = new Solution(&mesh, 0.0);
  
  // Initialize the weak formulation.
  WeakForm wf;
  info("Registering forms for the Newton's method.");
  wf.add_matrix_form(jac_form_vol, jac_form_vol_ord, HERMES_NONSYM, HERMES_ANY, sln_time_prev);
  wf.add_vector_form(res_form_vol, res_form_vol_ord, HERMES_ANY, sln_time_prev);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, space, is_linear);

  // Visualize the projection and mesh.
  ScalarView sview("Initial condition", new WinGeom(0, 0, 400, 350));
  sview.fix_scale_width(50);
  sview.show(sln_time_prev, HERMES_EPS_VERYHIGH);
  ScalarView eview("Temporal error", new WinGeom(405, 0, 400, 350));
  eview.fix_scale_width(50);
  eview.show(time_error_fn, HERMES_EPS_VERYHIGH);
  OrderView oview("Initial mesh", new WinGeom(810, 0, 350, 350));
  oview.show(space);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  info("Time step history will be saved to file time_step_history.dat.");

  // Time stepping loop:
  int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    bool is_linear = false;
    if (!rk_time_step(current_time, time_step, &bt, sln_time_prev, sln_time_new, time_error_fn, &dp, matrix_solver,
		      verbose, is_linear, NEWTON_TOL, NEWTON_MAX_ITER)) {
      info("Runge-Kutta time step failed, decreasing time step size from %g to %g days.", 
           time_step, time_step * time_step_dec);
           time_step *= time_step_dec;
           if (time_step < time_step_min) error("Time step became too small.");
	   continue;
    }

    // Show error function.
    char title[100];
    sprintf(title, "Temporal error, t = %g", current_time);
    eview.set_title(title);
    eview.show(time_error_fn, HERMES_EPS_VERYHIGH);

    // Calculate relative time stepping error and decide whether the 
    // time step can be accepted. If not, then the time step size is 
    // reduced and the entire time step repeated. If yes, then another
    // check is run, and if the relative error is very low, time step 
    // is increased.
    double rel_err_time = calc_norm(time_error_fn, HERMES_H1_NORM) / calc_norm(sln_time_new, HERMES_H1_NORM) * 100;
    info("rel_err_time = %g%%", rel_err_time);
    if (rel_err_time > time_tol_upper) {
      info("rel_err_time above upper limit %g%% -> decreasing time step from %g to %g days and repeating time step.", 
           time_tol_upper, time_step, time_step * time_step_dec);
      time_step *= time_step_dec;
      continue;
    }
    if (rel_err_time < time_tol_lower) {
      info("rel_err_time = below lower limit %g%% -> increasing time step from %g to %g days", 
           time_tol_lower, time_step, time_step * time_step_inc);
      time_step *= time_step_inc;
    }

    // Add entry to the timestep graph.
    time_step_graph.add_values(current_time, time_step);
    time_step_graph.save("time_step_history.dat");

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(sln_time_new, HERMES_EPS_VERYHIGH);
    oview.show(space);

    // Save complete Solution.
    char filename[100];
    sprintf(filename, "outputs/tsln_%f.dat", current_time);
    bool compress = false;   // Gzip compression not used as it only works on Linux.
    sln_time_new->save(filename, compress);
    info("Solution at time %g saved to file %s.", current_time, filename);

    // Save solution for the next time step.
    sln_time_prev->copy(sln_time_new);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  // Cleanup.
  delete space;
  delete sln_time_prev;
  delete sln_time_new;
  delete time_error_fn;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
