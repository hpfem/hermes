#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example solves adaptively a time-dependent coupled problem of heat and moisture 
// transfer in massive concrete walls of a nuclear reactor vessel (simplified axi-symmetric 
// geometry). 
//
// PDE: Lengthy. See the paper P. Solin, L. Dubcova, J. Kruis: Adaptive hp-FEM with Dynamical 
// Meshes for Transient Heat and Moisture Transfer Problems, J. Comput. Appl. Math. 233 (2010) 3103-3112.
//
// The following parameters can be changed:

const int P_INIT = 1;                             // Initial polynomial degrees.
const bool MULTI = true;                          // MULTI = true  ... use multi-mesh,
                                                  // MULTI = false ... use single-mesh.
                                                  // Note: In the single mesh option, the meshes are
                                                  // forced to be geometrically the same but the
                                                  // polynomial degrees can still vary.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_LEVEL = 1;                        // 1 = one layer of refinements is shaved off and poly degrees 
                                                  // of all elements reset to P_INIT; 2 = mesh reset to basemesh.  
                                                  // TODO: Add a third option where one layer will be taken off 
                                                  // and just one polynomial degree subtracted.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.3;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time step and simulation time.
const double TAU = 5.*24*60*60;                  // time step: 120 hours
const double SIMULATION_TIME = 100*TAU + 0.001;  // (seconds) physical time

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
const int BDY_SYMMETRY = 1;               
const int BDY_REACTOR_WALL = 2;           
const int BDY_EXTERIOR_WALL = 5;          

// Physical time in seconds.
double CURRENT_TIME = 0.0;

// Essential (Dirichlet) boundary condition values for T.
scalar essential_bc_values_T(double x, double y, double time)
{
  double current_reactor_temperature = TEMP_REACTOR_MAX;
  if (time < REACTOR_START_TIME) {
    current_reactor_temperature = TEMP_INITIAL +
      (time/REACTOR_START_TIME)*(TEMP_REACTOR_MAX - TEMP_INITIAL);
  }
  return current_reactor_temperature;
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

  // Enter boundary markers.
  BCTypes temp_bc_type, moist_bc_type;
  temp_bc_type.add_bc_dirichlet(BDY_REACTOR_WALL);
  temp_bc_type.add_bc_neumann(BDY_SYMMETRY);
  temp_bc_type.add_bc_newton(BDY_EXTERIOR_WALL);
  moist_bc_type.add_bc_neumann(Hermes::vector<int>(BDY_SYMMETRY, BDY_REACTOR_WALL));
  moist_bc_type.add_bc_newton(BDY_EXTERIOR_WALL);

  // Enter Dirichlet boundary values.
  BCValues bc_values(&CURRENT_TIME);
  bc_values.add_timedep_function(BDY_REACTOR_WALL, essential_bc_values_T);

  // Create H1 spaces with default shapesets.
  H1Space T_space(&T_mesh, &temp_bc_type, &bc_values, P_INIT);
  H1Space M_space(MULTI ? &M_mesh : &T_mesh, &moist_bc_type, P_INIT);

  // Define constant initial conditions.
  info("Setting initial conditions.");
  Solution T_prev_time(&T_mesh, TEMP_INITIAL);
  Solution M_prev_time(&M_mesh, MOIST_INITIAL);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_sym_0_1));
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_1_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_sym_1_0));
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, &T_prev_time);
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, &M_prev_time);
  wf.add_matrix_form_surf(0, 0, callback(bilinear_form_surf_0_0_ext), BDY_EXTERIOR_WALL);
  wf.add_matrix_form_surf(1, 1, callback(bilinear_form_surf_1_1_ext), BDY_EXTERIOR_WALL);
  wf.add_vector_form_surf(0, callback(linear_form_surf_0_ext), BDY_EXTERIOR_WALL);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1_ext), BDY_EXTERIOR_WALL);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions.
  Solution T_coarse, M_coarse, T_fine, M_fine;

  // Geometry and position of visualization windows.
  WinGeom* T_sln_win_geom = new WinGeom(0, 0, 300, 450);
  WinGeom* M_sln_win_geom = new WinGeom(310, 0, 300, 450);
  WinGeom* T_mesh_win_geom = new WinGeom(620, 0, 280, 450);
  WinGeom* M_mesh_win_geom = new WinGeom(910, 0, 280, 450);

  // Initialize views.
  ScalarView T_sln_view("Temperature", T_sln_win_geom);
  ScalarView M_sln_view("Moisture", M_sln_win_geom);
  OrderView T_order_view("Temperature mesh", T_mesh_win_geom);
  OrderView M_order_view("Moisture mesh", M_mesh_win_geom);

  // Show initial conditions.
  T_sln_view.show(&T_prev_time);
  M_sln_view.show(&M_prev_time);
  T_order_view.show(&T_space);
  M_order_view.show(&M_space);

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
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      if (UNREF_LEVEL == 1) {
        T_mesh.unrefine_all_elements();
        M_mesh.unrefine_all_elements();
      }
      else {
        T_mesh.copy(&basemesh);
        M_mesh.copy(&basemesh);
      }
      T_space.set_uniform_order(P_INIT);
      M_space.set_uniform_order(P_INIT);
    }

    // Spatial adaptivity loop. Note: T_prev_time and M_prev_time must not be changed during 
    // spatial adaptivity.
    int as = 1; 
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&T_space, &M_space));

      // Initialize matrix solver.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      bool is_linear = true;
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      dp->assemble(matrix, rhs);

      // Now we can deallocate the previous fine meshes.
      if(as > 1){ delete T_fine.get_mesh(); delete M_fine.get_mesh(); }

      // Solve the linear system of the reference problem. If successful, obtain the solutions.
      if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                              Hermes::vector<Solution *>(&T_fine, &M_fine));
      else error ("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&T_space, &M_space), 
                                   Hermes::vector<Solution *>(&T_fine, &M_fine), 
                                   Hermes::vector<Solution *>(&T_coarse, &M_coarse), matrix_solver); 

      // Registering custom forms for error calculation.
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&T_space, &M_space));
      adaptivity->set_error_form(0, 0, callback(bilinear_form_sym_0_0));
      adaptivity->set_error_form(0, 1, callback(bilinear_form_sym_0_1));
      adaptivity->set_error_form(1, 0, callback(bilinear_form_sym_1_0));
      adaptivity->set_error_form(1, 1, callback(bilinear_form_sym_1_1));

      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&T_coarse, &M_coarse), 
                                 Hermes::vector<Solution *>(&T_fine, &M_fine)) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space::get_num_dofs(Hermes::vector<Space *>(&T_space, &M_space)), 
                            Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) 
        done = true;
      else 
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(Hermes::vector<Space *>(&T_space, &M_space)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      delete ref_spaces;
      delete dp;
      
      // Increase counter.
      as++;
    }
    while (done == false);

    // Update time.
    CURRENT_TIME += TAU;

    // Show new coarse meshes and solutions.
    char title[100];
    sprintf(title, "Temperature, t = %g days", CURRENT_TIME/3600./24);
    T_sln_view.set_title(title);
    T_sln_view.show(&T_coarse);
    sprintf(title, "Moisture, t = %g days", CURRENT_TIME/3600./24);
    M_sln_view.set_title(title);
    M_sln_view.show(&M_coarse);
    T_order_view.show(&T_space);
    M_order_view.show(&M_space);

    // Save fine mesh solutions for the next time step.
    T_prev_time.copy(&T_fine);
    M_prev_time.copy(&M_fine);

    ts++;
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
