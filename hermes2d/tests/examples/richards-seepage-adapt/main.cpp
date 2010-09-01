#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example "richards-seepage-adapt" works correctly.

const int P_INIT = 1;                      // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 0;            // Number of initial mesh refinements towards the top edge.
const int TIME_INTEGRATION = 2;            // 1... implicit Euler, 2... Crank-Nicolson.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;   // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See the User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Newton's method
const double NEWTON_TOL_COARSE = 0.0001;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.0005;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 50;            // Maximum allowed number of Newton iterations.
const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method.

// Problem parameters.
const double TAU = 5e-3;                   // Time step.
const double STARTUP_TIME = 1.1e-2;        // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 5.0;                // Time interval length.
double TIME = 0;                           // Global time variable initialized with first time step.
double H_INIT = -9.5;                      // Initial pressure head.
double H_ELEVATION = 5.2;

double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
double K_S_4 = 1.061;

double ALPHA_1 = 0.01;
double ALPHA_3 = 0.005;
double ALPHA_2 = 0.01;
double ALPHA_4 = 0.05;

double THETA_R_1 = 0.1020;
double THETA_R_2 = 0.09849;
double THETA_R_3 = 0.08590;
double THETA_R_4 = 0.08590;

double THETA_S_1 = 0.4570;
double THETA_S_2 = 0.4510;
double THETA_S_3 = 0.4650;
double THETA_S_4 = 0.5650;

double N_1 = 1.982;
double N_2 = 1.632; 
double N_3 = 5.0;
double N_4 = 5.0;

double M_1 = 0.49546;
double M_2 = 0.38726;
double M_3 = 0.8;
double M_4 = 0.8;

double Q_MAX_VALUE = 0.07;         // Maximum value, used in function q_function(); 
double q_function() {
  if (STARTUP_TIME > TIME) return Q_MAX_VALUE * TIME / STARTUP_TIME;
  else return Q_MAX_VALUE;
}

double STORATIVITY = 0.05;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Material properties.
bool is_in_mat_1(double x, double y) {
  if (y >= -0.5) return true;
  else return false; 
}

bool is_in_mat_2(double x, double y) {
  if (y >= -1.0 && y < -0.5) return true;
  else return false; 
}

bool is_in_mat_4(double x, double y) {
  if (x >= 1.0 && x <= 3.0 && y >= -2.5 && y < -1.5) return true;
  else return false; 
}

bool is_in_mat_3(double x, double y) {
  if (!is_in_mat_1(x, y) && !is_in_mat_2(x, y) && !is_in_mat_4(x, y)) return true;
  else return false; 
}

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
int BDY_1 = 1;
int BDY_3 = 3;
int BDY_4 = 4;
int BDY_6 = 6;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == BDY_3) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Initial condition.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = 0;
  dy = -1;
  return -y + H_INIT;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  if (STARTUP_TIME > TIME) return -y + H_INIT + TIME/STARTUP_TIME*H_ELEVATION;
  else return -y + H_INIT + H_ELEVATION;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(3, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution sln, ref_sln, u_prev_time;

  // Adapt mesh to represent initial condition with given accuracy.
  int proj_norm = 1;  // H1 norm.
  bool verbose = false; 
  double err_stop_init_cond = 0.1 * ERR_STOP; 
  adapt_to_exact_function(&space, proj_norm, init_cond, &selector, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY, ERR_STOP, NDOF_STOP, 
                          verbose, &u_prev_time);

  // Assign initial condition to mesh.
  u_prev_time.set_exact(&mesh, init_cond);
  Vector *coeff_vec = new AVector(ndof);

  // Calculating initial vector for Newton.
  info("Projecting initial condition to obtain coefficient vector for Newton on coarse mesh.");
  project_global(&space, H2D_H1_NORM, &u_prev_time, &u_prev_time, coeff_vec);

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_form_vol_euler, jac_form_vol_ord, H2D_UNSYM, H2D_ANY, 
                       Tuple<MeshFunction*>(&u_prev_time));
    wf.add_matrix_form_surf(jac_form_surf_1_euler, jac_form_surf_1_ord, BDY_1);
    wf.add_matrix_form_surf(jac_form_surf_4_euler, jac_form_surf_4_ord, BDY_4);
    wf.add_matrix_form_surf(jac_form_surf_6_euler, jac_form_surf_6_ord, BDY_6);
    wf.add_vector_form(res_form_vol_euler, res_form_vol_ord, H2D_ANY, 
                       Tuple<MeshFunction*>(&u_prev_time));
    wf.add_vector_form_surf(res_form_surf_1_euler, res_form_surf_1_ord, BDY_1); 
    wf.add_vector_form_surf(res_form_surf_4_euler, res_form_surf_4_ord, BDY_4);
    wf.add_vector_form_surf(res_form_surf_6_euler, res_form_surf_6_ord, BDY_6);
  }
  else {
    wf.add_matrix_form(jac_form_vol_cranic, jac_form_vol_ord, H2D_UNSYM, H2D_ANY, 
                       Tuple<MeshFunction*>(&u_prev_time));
    wf.add_matrix_form_surf(jac_form_surf_1_cranic, jac_form_surf_1_ord, BDY_1);
    wf.add_matrix_form_surf(jac_form_surf_4_cranic, jac_form_surf_4_ord, BDY_4);
    wf.add_matrix_form_surf(jac_form_surf_6_cranic, jac_form_surf_6_ord, BDY_6); 
    wf.add_vector_form(res_form_vol_cranic, res_form_vol_ord, H2D_ANY, 
                       Tuple<MeshFunction*>(&u_prev_time));
    wf.add_vector_form_surf(res_form_surf_1_cranic, res_form_surf_1_ord, BDY_1, 
			    Tuple<MeshFunction*>( &u_prev_time));
    wf.add_vector_form_surf(res_form_surf_4_cranic, res_form_surf_4_ord, BDY_4, 
			    Tuple<MeshFunction*>(&u_prev_time));
    wf.add_vector_form_surf(res_form_surf_6_cranic, res_form_surf_6_ord, BDY_6, 
			    Tuple<MeshFunction*>(&u_prev_time));
  }
 
  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, MESH_REGULARITY);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }

    // Update the coefficient vector and u_prev_time.
    info("Projecting to obtain coefficient vector on coarse mesh.");
    project_global(&space, H2D_H1_NORM, &u_prev_time, &u_prev_time, coeff_vec);

    bool verbose = false;     // Print info during adaptivity.
    info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
    // The NULL pointers mean that we are not interested in visualization during the Newton's loop.
    solve_newton_adapt(&space, &wf, coeff_vec, matrix_solver, H2D_H1_NORM, &sln, &ref_sln,
                       Tuple<WinGeom *>(), Tuple<WinGeom *>(), &selector, &apt,
                       NEWTON_TOL_COARSE, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose);

    // Copy new time level reference solution into u_prev_time.
    u_prev_time.set_fe_solution(&space, coeff_vec);
  }

  // Waiting for test.
}
