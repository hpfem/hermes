#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 23-heat-and-moisture-adapt works correctly.

const int P_INIT = 1;                    // Initial polynomial degrees.
const bool MULTI = true;                 // MULTI = true  ... use multi-mesh,
                                         // MULTI = false ... use single-mesh.
                                         // Note: In the single mesh option, the meshes are
                                         // forced to be geometrically the same but the
                                         // polynomial degrees can still vary.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;            // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Time step and simulation time.
const double TAU = 5.*24*60*60;                 // time step: 120 hours
const double SIMULATION_TIME = TAU + 0.001;     // (seconds) physical time

// Equation parameters.
const double c_TT = 2.18e+6;
const double d_TT = 2.1;
const double d_Tw = 2.37e-2;
const double k_TT = 25;
const double c_ww = 24.9;
const double d_wT = 1.78e-10;
const double d_ww = 3.02e-8;
const double k_ww = 1.84e-7;

// Initial and boundary conditions.
const double TEMP_INITIAL = 293.0;           // (Kelvins)
const double MOIST_INITIAL = 0.9;            // (dimensionless)
const double TEMP_EXTERIOR = 293.0;          // (Kelvins)
const double MOIST_EXTERIOR = 0.55;          // (dimensionless)
const double TEMP_REACTOR_MAX = 550.0;       // (Kelvins)
const double REACTOR_START_TIME = 3600*24;   // (seconds) how long does the reactor
                                             // need to warm up linearly from TEMP_INITIAL
                                             // to TEMP_REACTOR_MAX
// Materials and boundary markers.
const int MARKER_SYMMETRY = 1;               
const int MARKER_REACTOR_WALL = 2;           
const int MARKER_EXTERIOR_WALL = 5;          

// Physical time in seconds.
double CURRENT_TIME = 0.0;

// Boundary condition types.
BCType temp_bc_type(int marker)
  { return (marker == MARKER_REACTOR_WALL) ? BC_ESSENTIAL : BC_NATURAL; }

BCType moist_bc_type(int marker)
  { return BC_NATURAL; }

// Essential (Dirichlet) boundary condition values for T.
scalar essential_bc_values_T(int ess_bdy_marker, double x, double y)
{
  if (ess_bdy_marker == MARKER_REACTOR_WALL)
  {
    double current_reactor_temperature = TEMP_REACTOR_MAX;
    if (CURRENT_TIME < REACTOR_START_TIME) {
      current_reactor_temperature = TEMP_INITIAL +
        (CURRENT_TIME/REACTOR_START_TIME)*(TEMP_REACTOR_MAX - TEMP_INITIAL);
    }
    return current_reactor_temperature;
  }
  else return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh, T_mesh, M_mesh;
  H2DReader mloader;
  mloader.load("domain2.mesh", &basemesh);

  // Create temperature and moisture meshes.
  // This also initializes the multimesh hp-FEM.
  T_mesh.copy(&basemesh);
  M_mesh.copy(&basemesh);

  // Create H1 spaces with default shapesets.
  H1Space T_space(&T_mesh, temp_bc_type, essential_bc_values_T, P_INIT);
  H1Space M_space(MULTI ? &M_mesh : &T_mesh, moist_bc_type, NULL, P_INIT);

  // Define constant initial conditions.
  info("Setting initial conditions.");
  Solution T_prev, M_prev;
  T_prev.set_const(&T_mesh, TEMP_INITIAL);
  M_prev.set_const(&M_mesh, MOIST_INITIAL);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_sym_0_1));
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_1_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_sym_1_0));
  wf.add_vector_form(0, callback(linear_form_0), H2D_ANY, &T_prev);
  wf.add_vector_form(1, callback(linear_form_1), H2D_ANY, &M_prev);
  wf.add_matrix_form_surf(0, 0, callback(bilinear_form_surf_0_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_matrix_form_surf(1, 1, callback(bilinear_form_surf_1_1_ext), MARKER_EXTERIOR_WALL);
  wf.add_vector_form_surf(0, callback(linear_form_surf_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1_ext), MARKER_EXTERIOR_WALL);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  double to_be_processed = 0;
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY,
                          MESH_REGULARITY, to_be_processed, H2D_TOTAL_ERROR_REL, H2D_ELEMENT_ERROR_REL);
  apt.set_error_form(0, 0, callback(bilinear_form_sym_0_0));
  apt.set_error_form(0, 1, callback(bilinear_form_sym_0_1));
  apt.set_error_form(1, 0, callback(bilinear_form_sym_1_0));
  apt.set_error_form(1, 1, callback(bilinear_form_sym_1_1));

  // Solutions.
  Solution T_coarse, M_coarse, T_fine, M_fine;

  // Time stepping loop:
  bool verbose = true;  // Print info during adaptivity.
  double comp_time = 0.0;
  static int ts = 1;
  while (CURRENT_TIME < SIMULATION_TIME)
  {
    info("Simulation time = %g s (%d h, %d d, %d y)",
        (CURRENT_TIME + TAU), (int) (CURRENT_TIME + TAU) / 3600,
        (int) (CURRENT_TIME + TAU) / (3600*24), (int) (CURRENT_TIME + TAU) / (3600*24*364));

    // Uniform mesh derefinement.
    if (ts > 1) {
      T_mesh.copy(&basemesh);
      M_mesh.copy(&basemesh);
      T_space.set_uniform_order(P_INIT);
      M_space.set_uniform_order(P_INIT);
    }

    // Adaptivity loop.
    solve_linear_adapt(Tuple<Space *>(&T_space, &M_space), &wf, NULL, matrix_solver,
                       Tuple<int>(H2D_H1_NORM, H2D_H1_NORM),
                       Tuple<Solution *>(&T_coarse, &M_coarse),
                       Tuple<Solution *>(&T_fine, &M_fine),
                       Tuple<WinGeom *>(), Tuple<WinGeom *>(),// Do not show solutions or meshes.
                       Tuple<RefinementSelectors::Selector *> (&selector, &selector), &apt,
                       verbose);

    // Update time.
    CURRENT_TIME += TAU;

    // Save fine mesh solutions for the next time step.
    T_prev.copy(&T_fine);
    M_prev.copy(&M_fine);

    ts++;
  }
  info("Coordinate (  5.0,  4.0) T value = %lf", T_prev.get_pt_value(5.0, 4.0));
  info("Coordinate ( 12.0,  4.0) T value = %lf", T_prev.get_pt_value(12.0, 4.0));
  info("Coordinate ( 12.0, 17.0) T value = %lf", T_prev.get_pt_value(12.0, 17.0));
  info("Coordinate ( 12.0, 29.0) T value = %lf", T_prev.get_pt_value(12.0, 29.0));
  info("Coordinate (  5.0, 29.0) T value = %lf", T_prev.get_pt_value(5.0, 29.0));

  info("Coordinate (  5.0,  4.0) M value = %lf", M_prev.get_pt_value(5.0, 4.0));
  info("Coordinate ( 12.0,  4.0) M value = %lf", M_prev.get_pt_value(12.0, 4.0));
  info("Coordinate ( 12.0, 17.0) M value = %lf", M_prev.get_pt_value(12.0, 17.0));
  info("Coordinate ( 12.0, 29.0) M value = %lf", M_prev.get_pt_value(12.0, 29.0));
  info("Coordinate (  5.0, 29.0) M value = %lf", M_prev.get_pt_value(5.0, 29.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(T_prev.get_pt_value(5.0, 4.0) - 294.128000) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 4.0) - 293.091307) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 17.0) - 297.046535) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 29.0) - 293.666772) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(5.0, 29.0) - 308.780420) > eps) {
    success = 0;
  }

  if (fabs(M_prev.get_pt_value(5.0, 4.0) - 0.900089) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(12.0, 4.0) - 0.899224) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(12.0, 17.0) - 0.899553) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(12.0, 29.0) - 0.899292) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(5.0, 29.0) - 0.900779) > eps) {
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
