#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 24-heat-and-moisture-adapt works correctly.

const int P_INIT = 1;                             // Initial polynomial degrees.
const bool MULTI = true;                          // MULTI = true  ... use multi-mesh,
                                                  // MULTI = false ... use single-mesh.
                                                  // Note: In the single mesh option, the meshes are
                                                  // forced to be geometrically the same but the
                                                  // polynomial degrees can still vary.
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
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time step and simulation time.
const double TAU = 5.*24*60*60;                   // time step: 120 hours
const double SIMULATION_TIME = TAU + 0.001;       // (seconds) physical time

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
const double TEMP_INITIAL = 293.0;                // (Kelvins)
const double MOIST_INITIAL = 0.9;                 // (dimensionless)
const double TEMP_EXTERIOR = 293.0;               // (Kelvins)
const double MOIST_EXTERIOR = 0.55;               // (dimensionless)
const double TEMP_REACTOR_MAX = 550.0;            // (Kelvins)
const double REACTOR_START_TIME = 3600*24;        // (seconds) how long does the reactor
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
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

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
  Solution T_prev, M_prev;
  T_prev.set_const(&T_mesh, TEMP_INITIAL);
  M_prev.set_const(&M_mesh, MOIST_INITIAL);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_sym_0_1));
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_1_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_sym_1_0));
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, &T_prev);
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, &M_prev);
  wf.add_matrix_form_surf(0, 0, callback(bilinear_form_surf_0_0_ext), BDY_EXTERIOR_WALL);
  wf.add_matrix_form_surf(1, 1, callback(bilinear_form_surf_1_1_ext), BDY_EXTERIOR_WALL);
  wf.add_vector_form_surf(0, callback(linear_form_surf_0_ext), BDY_EXTERIOR_WALL);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1_ext), BDY_EXTERIOR_WALL);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

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
      info("Global mesh derefinement.");
      T_mesh.copy(&basemesh);
      M_mesh.copy(&basemesh);
      T_space.set_uniform_order(P_INIT);
      M_space.set_uniform_order(P_INIT);
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&T_space, &M_space));

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      bool is_linear = true;
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
      dp->assemble(matrix, rhs);

      // Now we can deallocate the previous fine meshes.
      if(as > 1){ delete T_fine.get_mesh(); delete M_fine.get_mesh(); }
      
      // Time measurement.
      cpu_time.tick();

      // Solve the linear system of the reference problem. If successful, obtain the solutions.
      if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                              Hermes::vector<Solution *>(&T_fine, &M_fine));
      else error ("Matrix solver failed.\n");
    
      // Time measurement.
      cpu_time.tick();

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&T_space, &M_space), Hermes::vector<Solution *>(&T_fine, &M_fine), 
                     Hermes::vector<Solution *>(&T_coarse, &M_coarse), matrix_solver); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&T_space, &M_space));
      adaptivity->set_error_form(0, 0, callback(bilinear_form_sym_0_0));
      adaptivity->set_error_form(0, 1, callback(bilinear_form_sym_0_1));
      adaptivity->set_error_form(1, 0, callback(bilinear_form_sym_1_0));
      adaptivity->set_error_form(1, 1, callback(bilinear_form_sym_1_1));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&T_coarse, &M_coarse), 
                                 Hermes::vector<Solution *>(&T_fine, &M_fine)) * 100;

      // Time measurement.
      cpu_time.tick();

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space::get_num_dofs(Hermes::vector<Space *>(&T_space, &M_space)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);
      
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

  int success = 1;
  double eps = 1e-5;
  if (fabs(T_prev.get_pt_value(5.0, 4.0) - 294.127903) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 4.0) - 293.091431) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 17.0) - 297.046520) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(12.0, 29.0) - 293.667172) > eps) {
    success = 0;
  }
  if (fabs(T_prev.get_pt_value(5.0, 29.0) - 308.780420) > eps) {
    success = 0;
  }

  if (fabs(M_prev.get_pt_value(5.0, 4.0) - 0.900090) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(12.0, 4.0) - 0.899229) > eps) {
    success = 0;
  }
  if (fabs(M_prev.get_pt_value(12.0, 17.0) - 0.898690) > eps) {
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
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
