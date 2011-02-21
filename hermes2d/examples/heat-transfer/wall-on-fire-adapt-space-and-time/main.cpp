#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

using namespace RefinementSelectors;

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

const int P_INIT = 3;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
double time_step = 20;                            // Time step in seconds.

// Spatial adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_LEVEL = 1;                        // 1 = one layer of refinements is shaved off and poly degrees
                                                  // of all elements reset to P_INIT; 2 = mesh reset to basemesh.  
                                                  // TODO: Add a third option where one layer will be taken off 
                                                  // and just one polynomial degree subtracted.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
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
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double SPACE_ERR_TOL = 1.0;                 // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Temporal adaptivity.
bool ADAPTIVE_TIME_STEP_ON = true;                // This flag decides whether adaptive time stepping will be done.
                                                  // The methods for the adaptive and fixed-step versions are set
                                                  // below. An embedded method must be used with adaptive time stepping. 
const double TIME_ERR_TOL_UPPER = 1.0;                // If rel. temporal error is greater than this threshold, decrease time 
                                                  // step size and repeat time step.
const double TIME_ERR_TOL_LOWER = 0.5;                // If rel. temporal error is less than this threshold, increase time step
                                                  // but do not repeat time step (this might need further research).
const double TIME_STEP_INC_RATIO = 1.1;           // Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_DEC_RATIO = 0.8;           // Time step decrease ratio (applied when rel. temporal error is too large).

// Newton's method.
const double NEWTON_TOL_COARSE = 0.001;           // Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.005;             // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

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
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

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
const double HEATCAP = 1020;       // Heat capacity.
const double RHO = 2200;           // Material density.
const double T_FINAL = 4000;       // Length of time interval in seconds.

// Problem-specific functions.
#include "extras.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable* bt = new ButcherTable(butcher_table_type);
  if (bt->is_explicit()) info("Using a %d-stage explicit R-K method.", bt->get_size());
  if (bt->is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt->get_size());
  if (bt->is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt->get_size());

  // Turn off adaptive time stepping if R-K method is not embedded.
  if (bt->is_embedded() == false && ADAPTIVE_TIME_STEP_ON == true) {
    warn("R-K method not embedded, turning off adaptive time stepping.");
    ADAPTIVE_TIME_STEP_ON = false;
  }

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("wall.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(BDY_BOTTOM, INIT_REF_NUM_BDY);
  mesh.copy(&basemesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_RIGHT, BDY_LEFT));
  bc_types.add_bc_newton(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP));

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, NULL, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Convert initial condition into a Solution.
  Solution sln_prev_time(&mesh, TEMP_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(stac_jacobian_vol, stac_jacobian_vol_ord, HERMES_NONSYM, HERMES_ANY, &sln_prev_time);
  wf.add_vector_form(stac_residual_vol, stac_residual_vol_ord, HERMES_ANY, &sln_prev_time);
  wf.add_matrix_form_surf(stac_jacobian_bottom, stac_jacobian_bottom_ord, BDY_BOTTOM, &sln_prev_time);
  wf.add_vector_form_surf(stac_residual_bottom, stac_residual_bottom_ord, BDY_BOTTOM, &sln_prev_time);
  wf.add_matrix_form_surf(stac_jacobian_top, stac_jacobian_top_ord, BDY_TOP, &sln_prev_time);
  wf.add_vector_form_surf(stac_residual_top, stac_residual_top_ord, BDY_TOP, &sln_prev_time);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  char title[100];
  ScalarView sln_view("Initial condition", new WinGeom(0, 0, 1500, 360));
  OrderView ordview("Initial mesh", new WinGeom(0, 410, 1500, 360));
  ScalarView time_error_view("Temporal error", new WinGeom(0, 800, 1500, 360));
  time_error_view.fix_scale_width(40);
  ScalarView space_error_view("Spatial error", new WinGeom(0, 1220, 1500, 360));
  space_error_view.fix_scale_width(40);
  sln_view.show(&sln_prev_time, HERMES_EPS_VERYHIGH);
  ordview.show(&space);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  if (ADAPTIVE_TIME_STEP_ON) info("Time step history will be saved to file time_step_history.dat.");

  // Time stepping loop:
  double current_time = 0; int ts = 1;
  do 
  {
    info("Begin time step %d.", ts);
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      if (UNREF_LEVEL == 1) mesh.unrefine_all_elements();
      else mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      ndof = Space::get_num_dofs(&space);
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be 
    // changed during spatial adaptivity. 
    Solution ref_sln;
    Solution time_error_fn(&mesh);
    bool done = false; int as = 1;
    double err_est;
    do {
      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = construct_refined_space(&space);

      OGProjection::project_global(ref_space, Hermes::vector<Solution *>(&sln_prev_time), 
                     Hermes::vector<Solution *>(&sln_prev_time), matrix_solver);
      
      delete ref_sln.get_mesh();
      
      // Initialize discrete problem on reference mesh.
      DiscreteProblem* ref_dp = new DiscreteProblem(&wf, ref_space);

      // Runge-Kutta step on the fine mesh.
      info("Runge-Kutta time step on fine mesh (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt->get_size());
      bool verbose = true;
      bool is_linear = false;
      if (!rk_time_step(current_time, time_step, bt, &sln_prev_time, &ref_sln, bt->is_embedded() ? &time_error_fn : NULL,
                        ref_dp, matrix_solver, verbose, is_linear, NEWTON_TOL_FINE, NEWTON_MAX_ITER)) {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      /* If ADAPTIVE_TIME_STEP_ON == true, estimate temporal error. 
         If too large or too small, then adjust it and restart the time step. */

      double rel_err_time;
      if (bt->is_embedded() == true) {
        info("Calculating temporal error estimate.");

        // Show temporal error.
        char title[100];
        sprintf(title, "Temporal error est, spatial adaptivity step %d", as);     
        time_error_view.set_title(title);
        time_error_view.show_mesh(false);
        time_error_view.show(&time_error_fn, HERMES_EPS_VERYHIGH);

        rel_err_time = calc_norm(&time_error_fn, HERMES_H1_NORM) / calc_norm(&ref_sln, HERMES_H1_NORM) * 100;
        if (ADAPTIVE_TIME_STEP_ON == false) info("rel_err_time: %g%%", rel_err_time);
      }

      if (ADAPTIVE_TIME_STEP_ON) {
        if (rel_err_time > TIME_ERR_TOL_UPPER) {
          info("rel_err_time %g%% is above upper limit %g%%", rel_err_time, TIME_ERR_TOL_UPPER);
          info("Decreasing tau from %g to %g s and restarting time step.", 
               time_step, time_step * TIME_STEP_DEC_RATIO);
          time_step *= TIME_STEP_DEC_RATIO;
          delete ref_space;
          delete ref_dp;
          continue;
        }
        else if (rel_err_time < TIME_ERR_TOL_LOWER) {
          info("rel_err_time = %g%% is below lower limit %g%%", rel_err_time, TIME_ERR_TOL_UPPER);
          info("Increasing tau from %g to %g s.", time_step, time_step * TIME_STEP_INC_RATIO);
          time_step *= TIME_STEP_INC_RATIO;
        }
        else {
          info("rel_err_time = %g%% is in acceptable interval (%g%%, %g%%)", 
            rel_err_time, TIME_ERR_TOL_LOWER, TIME_ERR_TOL_UPPER);
        }

        // Add entry to time step history graph.
        time_step_graph.add_values(current_time, time_step);
        time_step_graph.save("time_step_history.dat");
      }

      /* Estimate spatial errors and perform mesh refinement */

      info("Spatial adaptivity step %d.", as);

      // Project the fine mesh solution onto the coarse mesh.
      Solution sln;
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Show spatial error.
      sprintf(title, "Spatial error est, spatial adaptivity step %d", as);  
      DiffFilter space_error_fn(Hermes::vector<MeshFunction*>(&ref_sln, &sln));   
      space_error_view.set_title(title);
      space_error_view.show_mesh(false);
      AbsFilter abs_sef(&space_error_fn);
      space_error_view.show(&abs_sef, HERMES_EPS_VERYHIGH);

      // Calculate element errors and spatial error estimate.
      info("Calculating spatial error estimate.");
      Adapt* adaptivity = new Adapt(&space);
      double err_rel_space = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_rel_space: %g%%", 
           Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_rel_space);

      // If err_est too large, adapt the mesh.
      if (err_rel_space < SPACE_ERR_TOL) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(&space) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
      
      // Clean up.
      delete adaptivity;
      delete ref_space;
      delete ref_dp;
    }
    while (done == false);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g s", current_time);
    sln_view.set_title(title);
    sln_view.show_mesh(false);
    sln_view.show(&ref_sln, HERMES_EPS_VERYHIGH);
    sprintf(title, "Mesh, time %g s", current_time);
    ordview.set_title(title);
    ordview.show(&space);

    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Clean up.
  delete bt;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
