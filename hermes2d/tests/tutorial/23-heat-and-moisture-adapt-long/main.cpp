#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of test example "23-heat-and-moisture-adapt-long": 
// function solve_linear_adapt() is not used.
// This test makes sure that example 23-heat-and-moisture-adapt-long works correctly.

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
const double SIMULATION_TIME = TAU + 0.001;    // (seconds) physical time

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
scalar T_essential_bc_values(int ess_bdy_marker, double x, double y)
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
  H1Space T_space(&T_mesh, temp_bc_type, T_essential_bc_values, P_INIT);
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

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err, graph_time_dof;

  // Solutions.
  Solution T_coarse, M_coarse, T_fine, M_fine;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, get_num_dofs(Tuple<Space *>(&T_space, &M_space)), mat, rhs, solver);

  // Time stepping loop:
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

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);
      info("Solving on reference mesh.");

      // Construct globally refined reference mesh
      // and setup reference space.
      int order_increase = 1;

      Mesh *T_ref_mesh = new Mesh();
      T_ref_mesh->copy(T_space.get_mesh());
      T_ref_mesh->refine_all_elements();
      Space* T_ref_space = T_space.dup(T_ref_mesh);
      T_ref_space->copy_orders(&T_space, order_increase);

      Mesh *M_ref_mesh = new Mesh();
      M_ref_mesh->copy(M_space.get_mesh());
      M_ref_mesh->refine_all_elements();
      Space* M_ref_space = M_space.dup(M_ref_mesh);
      M_ref_space->copy_orders(&M_space, order_increase);

      // Solve the reference problem.
      solve_linear(Tuple<Space *>(T_ref_space, M_ref_space), &wf, matrix_solver,
                   Tuple<Solution *>(&T_fine, &M_fine));

      // Project the reference solution on the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      // NULL means that we do not want to know the resulting coefficient vector.
      project_global(Tuple<Space *>(&T_space, &M_space),
                     Tuple<int>(H2D_H1_NORM, H2D_H1_NORM),
                     Tuple<MeshFunction *>(&T_fine, &M_fine),
                     Tuple<Solution *>(&T_coarse, &M_coarse), NULL);

      // Time measurement.
      cpu_time.tick();

      // Skip visualization time.
      cpu_time.tick(HERMES_SKIP);

      // Initialize the adaptivity module, set the coarse and fine mesh 
      // solutions, and set the error form.
      Adapt hp(Tuple<Space *>(&T_space, &M_space),
               Tuple<int>(H2D_H1_NORM, H2D_H1_NORM));
      hp.set_solutions(Tuple<Solution *>(&T_coarse, &M_coarse),
                       Tuple<Solution *>(&T_fine, &M_fine));
      hp.set_error_form(0, 0, callback(bilinear_form_sym_0_0));
      hp.set_error_form(0, 1, callback(bilinear_form_sym_0_1));
      hp.set_error_form(1, 0, callback(bilinear_form_sym_1_0));
      hp.set_error_form(1, 1, callback(bilinear_form_sym_1_1));

      // Calculate element errors.
      info("Calculating error (est).");
      hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL);

      // Calculate error estimate for each solution component.
      double T_err_est_abs = calc_abs_error(&T_coarse, &T_fine, H2D_H1_NORM);
      double T_norm_est = calc_norm(&T_fine, H2D_H1_NORM);
      double M_err_est_abs = calc_abs_error(&M_coarse, &M_fine, H2D_H1_NORM);
      double M_norm_est = calc_norm(&M_fine, H2D_H1_NORM);
      double err_est_abs_total = sqrt(T_err_est_abs*T_err_est_abs + M_err_est_abs*M_err_est_abs);
      double norm_est_total = sqrt(T_norm_est*T_norm_est + M_norm_est*M_norm_est);
      double err_est_rel_total = err_est_abs_total / norm_est_total * 100.;
      space_err_est = err_est_rel_total;

      // Report results.
      info("ndof[0]: %d, ref_ndof[0]: %d, err_est_rel[0]: %g%%",
           get_num_dofs(&T_space), get_num_dofs(T_ref_space),
           T_err_est_abs/T_norm_est*100);
      info("ndof[1]: %d, ref_ndof[1]: %d, err_est_rel[1]: %g%%",
           get_num_dofs(&M_space), get_num_dofs(M_ref_space),
           M_err_est_abs/M_norm_est*100);
      info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%",
           get_num_dofs(Tuple<Space *>(&T_space, &M_space)),
           get_num_dofs(Tuple<Space *>(T_ref_space, M_ref_space)), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (space_err_est > ERR_STOP) {
        info("Adapting coarse meshes.");
        done = hp.adapt(Tuple<RefinementSelectors::Selector *>(&selector, &selector), 
                        THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (get_num_dofs(Tuple<Space *>(&T_space, &M_space)) >= NDOF_STOP) done = true;
      }
      else done = true;

      as++;
    }
    while (!done);

    // Add entries to convergence graphs.
    graph_time_err.add_values(ts*TAU, space_err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(ts*TAU, get_num_dofs(Tuple<Space *>(&T_space, &M_space)));
    graph_time_dof.save("time_dof.dat");

    // Update time.
    CURRENT_TIME += TAU;

    // Save solutions for the next time step.
    T_prev.copy(&T_fine);
    M_prev.copy(&M_fine);

    ts++;
  }

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  printf("Success!\n");
  return ERROR_SUCCESS;

}
