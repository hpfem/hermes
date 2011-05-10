#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example is derived from example P03-timedep/03-newton-heat-ie
//  and it shows how automatic adaptivity in space can be combined with 
//  the implicit Euler method for a nonlinear time-dependent PDE. The
//  example uses fixed time step size.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity:
//
//  du/dt - div[lambda(u)grad u] - f = 0.
//
//  Domain: square (-10,10)^2.
//
//  BC:  Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double time_step = 0.5;                     // Time step. 
const double T_FINAL = 2.0;                       // Time interval length.

// Adaptivity
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_METHOD = 3;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
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
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;                   // Maximum allowed number of Newton iterations.

const double ALPHA = 4.0;                         // For the nonlinear thermal conductivity.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements(0, true);
  mesh.copy(&basemesh);
  
  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;

  // Convert initial condition into a Solution.
  InitialSolutionHeatTransfer sln_prev_time(&mesh);

  // Initialize the weak formulation
  WeakFormHeatTransferNewtonTimedep wf(ALPHA, time_step, &sln_prev_time);

  // Initialize the discrete problem.
  DiscreteProblem dp_coarse(&wf, &space);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(445, 0, 410, 350));
  view.show(&sln_prev_time);
  ordview.show(&space);
  
  // Time stepping loop.
  double current_time = time_step; int ts = 1;
  do 
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                //space.adjust_element_order(-1, P_INIT);
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof = Space::get_num_dofs(&space);
    }

    // The following is done only in the first time step, 
    // when the nonlinear problem was never solved before.
    if (ts == 1) {
      // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
      SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
      Vector* rhs_coarse = create_vector(matrix_solver);
      Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);
      scalar* coeff_vec_coarse = new scalar[ndof];

      // Calculate initial coefficient vector for Newton on the coarse mesh.
      info("Projecting initial condition to obtain coefficient vector on coarse mesh.");
      OGProjection::project_global(&space, &sln_prev_time, coeff_vec_coarse, matrix_solver);

      // Newton's loop on the coarse mesh.
      info("Solving on coarse mesh:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes2d.solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
	  jacobian_changed, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");
      Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

      // Cleanup after the Newton loop on the coarse mesh.
      delete matrix_coarse;
      delete rhs_coarse;
      delete solver_coarse;
      delete [] coeff_vec_coarse;
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be changed during spatial adaptivity. 
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = Space::construct_refined_space(&space);

      // Initialize matrix solver.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
      scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];

      // Initialize discrete problem on reference mesh.
      DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (ts == 1 && as == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on fine mesh.");
        OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
      }
      else {
        info("Projecting last fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
      }

      // Now we can deallocate the previous fine mesh.
      if(as > 1) delete ref_sln.get_mesh();

      // Newton's loop on the fine mesh.
      info("Solving on fine mesh:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes2d.solve_newton(coeff_vec, dp, solver, matrix, rhs, 
	  jacobian_changed, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Store the result in ref_sln.
      Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

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
      
      // Visualize the solution and mesh.
      char title[100];
      sprintf(title, "Solution, time %g", current_time);
      view.set_title(title);
      view.show_mesh(false);
      view.show(&ref_sln);
      sprintf(title, "Mesh, time %g", current_time);
      ordview.set_title(title);
      ordview.show(&space);

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      delete ref_space;
      delete dp;
      delete [] coeff_vec;
    }
    while (done == false);

    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
