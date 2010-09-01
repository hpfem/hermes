#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example "ns-timedep-adapt" works correctly.

const bool SOLVE_ON_COARSE_MESH = false; // true... Newton is done on coarse mesh in every adaptivity step.
                                         // false...Newton is done on coarse mesh only once, then projection
                                         // of the fine mesh solution to coarse mesh is used.
const int INIT_REF_NUM = 0;              // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;          // Number of initial mesh refinements towards boundary.
#define PRESSURE_IN_L2                   // If this is defined, the pressure is approximated using
                                         // discontinuous L2 elements (making the velocity discreetely
                                         // divergence-free, more accurate than using a continuous
                                         // pressure approximation). Otherwise the standard continuous
                                         // elements are used. The results are striking - check the
                                         // tutorial for comparisons.
const int P_INIT_VEL = 2;                // Initial polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;           // Initial polynomial degree for pressure
                                         // Note: P_INIT_VEL should always be greater than
                                         // P_INIT_PRESSURE because of the inf-sup condition

// Adaptivity
const int UNREF_FREQ = 1;        // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;    // This is a quantitative parameter of the adapt(...) function and
                                 // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          // Adaptive strategy:
                                 // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                 //   error is processed. If more elements have similar errors, refine
                                 //   all to keep the mesh symmetric.
                                 // STRATEGY = 1 ... refine all elements whose error is larger
                                 //   than THRESHOLD times maximum element error.
                                 // STRATEGY = 2 ... refine all elements whose error is larger
                                 //   than THRESHOLD.
                                 // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;  // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Used Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters
const double RE = 200.0;             // Reynolds number.
const double VEL_INLET = 1.0;        // Inlet velocity (reached after STARTUP_TIME).
const double STARTUP_TIME = 1.0;     // During this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant.
const double TAU = 0.01;             // Time step.
const double T_FINAL = 30000.0;      // Time interval length.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

// Geometry
const double H = 5;                  // Domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

// Boundary markers.
int bdy_bottom = 1;
int bdy_right  = 2;
int bdy_top = 3;
int bdy_left = 4;
int bdy_obstacle = 5;

// Current time (defined as global since needed in weak forms)
double TIME = 0;

// Boundary condition types for x-velocity
BCType xvel_bc_type(int marker) {
  if (marker == bdy_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

// Boundary condition values for x-velocity
scalar essential_bc_values_xvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_left) {
    // time-dependent inlet velocity (parabolic profile)
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

// Essential (Dirichlet) boundary condition values for y-velocity.
scalar essential_bc_values_yvel(int ess_bdy_marker, double x, double y) 
{
  return 0;
}

// Boundary condition types for y-velocity
BCType yvel_bc_type(int marker) {
  if (marker == bdy_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

BCType p_bc_type(int marker)
  { return BC_NONE; }

// Weak forms
#include "forms.cpp"

void mag(int n, scalar* a, scalar* dadx, scalar* dady,
                scalar* b, scalar* dbdx, scalar* dbdy,
                scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = sqrt(sqr(a[i]) + sqr(b[i]));
    outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dadx[i] + 2.0 * b[i] * dbdx[i]);
    outdy[i] = (0.5 / out[i]) * (2.0 * a[i] * dady[i] + 2.0 * b[i] * dbdy[i]);
  }
}

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh basemesh, mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);  // Master mesh.

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(bdy_obstacle, INIT_REF_NUM_BDY, false); // 'true' stands for anisotropic refinements,
  basemesh.refine_towards_boundary(bdy_top, INIT_REF_NUM_BDY, true);       // 'false' for isotropic.
  basemesh.refine_towards_boundary(bdy_bottom, INIT_REF_NUM_BDY, true);
  mesh.copy(&basemesh);

  // Create spaces with default shapesets. 
  H1Space* xvel_space = new H1Space(&mesh, xvel_bc_type, essential_bc_values_xvel, P_INIT_VEL);
  H1Space* yvel_space = new H1Space(&mesh, yvel_bc_type, essential_bc_values_yvel, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space* p_space = new L2Space(&mesh, P_INIT_PRESSURE);
#else
  H1Space* p_space = new H1Space(&mesh, p_bc_type, NULL, P_INIT_PRESSURE);
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = get_num_dofs(Tuple<Space *>(xvel_space, yvel_space, p_space));
  info("ndof = %d.", ndof);

  // Define projection norms.
  int vel_proj_norm = H2D_H1_NORM;
#ifdef PRESSURE_IN_L2
  int p_proj_norm = H2D_L2_NORM;
#else
  int p_proj_norm = H2D_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
//  Solution xvel_fine, yvel_fine, p_fine;
  Solution xvel_sln, yvel_sln, p_sln;
  Solution xvel_ref_sln, yvel_ref_sln, p_ref_sln;
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;

  // Define initial conditions on the coarse mesh.
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
  wf.add_matrix_form(0, 0, callback(newton_bilinear_form_unsym_0_0), H2D_UNSYM, H2D_ANY);
  wf.add_matrix_form(0, 1, callback(newton_bilinear_form_unsym_0_1), H2D_UNSYM, H2D_ANY);
  wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
  wf.add_matrix_form(1, 0, callback(newton_bilinear_form_unsym_1_0), H2D_UNSYM, H2D_ANY);
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
  wf.add_matrix_form(1, 1, callback(newton_bilinear_form_unsym_1_1), H2D_UNSYM, H2D_ANY);
  wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
  wf.add_vector_form(0, callback(newton_F_0), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
  wf.add_vector_form(1, callback(newton_F_1), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
  wf.add_vector_form(2, callback(newton_F_2), H2D_ANY);

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, MESH_REGULARITY);
  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  // Assign initial condition to mesh.
  Vector *coeff_vec = new AVector(ndof);

  // Time-stepping loop:
  char title[100];
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      xvel_space->set_uniform_order(P_INIT_VEL);
      yvel_space->set_uniform_order(P_INIT_VEL);
      p_space->set_uniform_order(P_INIT_PRESSURE);
    }

    // Update the coefficient vector and u_prev_time.
    info("Projecting to obtain coefficient vector on coarse mesh.");
    project_global(Tuple<Space *>(xvel_space, yvel_space, p_space),
                   Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm),
                   Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                   Tuple<Solution*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                   coeff_vec);

    // Adaptivity loop (in space):
    bool verbose = true;     // Print info during adaptivity.
    info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
    // The NULL pointers mean that we are not interested in visualization during the Newton's loop.
    solve_newton_adapt(Tuple<Space *>(xvel_space, yvel_space, p_space), &wf, coeff_vec, matrix_solver, 
                       Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm), 
                       Tuple<Solution *>(&xvel_sln, &yvel_sln, &p_sln),
                       Tuple<Solution *>(&xvel_ref_sln, &yvel_ref_sln, &p_ref_sln),
                       Tuple<WinGeom *>(), Tuple<WinGeom *>(), 
                       Tuple<RefinementSelectors::Selector *>(&selector, &selector, &selector), &apt,
                       NEWTON_TOL_COARSE, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose);

    // Copy new time level reference solution into prev_time.
    xvel_prev_time.set_fe_solution(xvel_space, coeff_vec);
    yvel_prev_time.set_fe_solution(yvel_space, coeff_vec);
    p_prev_time.set_fe_solution(p_space, coeff_vec);
  }

  ndof = get_num_dofs(Tuple<Space *>(xvel_space, yvel_space, p_space));
  info("ndof = %d", ndof);

  // Waiting for tests.
}
