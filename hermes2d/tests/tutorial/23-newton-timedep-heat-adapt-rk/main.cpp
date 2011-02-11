#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 23-newton-timedep-heat-adapt-rk works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double time_step = 0.5;                     // Time step. 
const double T_FINAL = time_step*2;               // Time interval length.

// Adaptivity
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Newton's method
const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;                   // Maximum allowed number of Newton iterations.

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
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/100.;
  dy = (x+10)/100.;
  return (x+10)*(y+10)/100.;
}

// Initial condition.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Convert initial condition into a Solution.
  Solution* sln_prev_time = new Solution(&mesh, init_cond);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  sln.copy(sln_prev_time);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian), HERMES_NONSYM, HERMES_ANY, sln_prev_time);
  wf.add_vector_form(callback(stac_residual), HERMES_ANY, sln_prev_time);

  // Initialize the discrete problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, &space, is_linear);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  
  // Time stepping loop.
  double current_time = time_step; int ts = 1;
  do 
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      if (UNREF_LEVEL == 1) mesh.unrefine_all_elements();
      else mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      ndof = Space::get_num_dofs(&space);
    }

    // The following is done only in the first time step, 
    // when the nonlinear problem was never solved before.
    if (ts == 1) {
      // Runge-Kutta step on the coarse mesh.
      bool verbose = true;
      bool is_linear = false;
      info("Runge-Kutta time step on coarse mesh (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
      if (!rk_time_step(current_time, time_step, &bt, &sln, &space, &dp_coarse, matrix_solver,
		        verbose, is_linear, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be changed during spatial adaptivity. 
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = construct_refined_space(&space);

      // Initialize discrete problem on reference mesh.
      DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);

      // Calculate initial solution on the fine mesh.
      if (ts == 1 && as == 1) {
        info("Projecting coarse mesh solution to obtain initial solution on fine mesh.");
        OGProjection::project_global(ref_space, &sln, &ref_sln, matrix_solver);
      }
      else {
        info("Projecting last fine mesh solution to obtain initial solution on new fine mesh.");
        OGProjection::project_global(ref_space, &ref_sln, &ref_sln, matrix_solver);
      }

      // Now we can deallocate the previous fine mesh.
      //if(as > 1) delete ref_sln.get_mesh();

      // Runge-Kutta step on the fine mesh.
      info("Runge-Kutta time step on fine mesh (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
      bool verbose = true;
      bool is_linear = false;
      if (!rk_time_step(current_time, time_step, &bt, &ref_sln, ref_space, dp, matrix_solver,
		        verbose, is_linear, NEWTON_TOL_FINE, NEWTON_MAX_ITER)) {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(&space);
      double err_est_rel_total = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
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
      delete dp;
    }
    while (done == false);

    // Copy last reference solution into sln_prev_time.
    sln_prev_time->copy(&ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Clean up.
  delete sln_prev_time;

  ndof = Space::get_num_dofs(&space);

  printf("ndof allowed = %d\n", 160);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 160) {      // ndofs was 153 at the time this test was created.
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
