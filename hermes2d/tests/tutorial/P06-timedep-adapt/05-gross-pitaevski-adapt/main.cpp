#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 25-newton-timedep-gp-adapt works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform refinements.
const int P_INIT = 2;                             // Initial polynomial degree.
const int TIME_DISCR = 2;                         // 1 for implicit Euler, 2 for Crank-Nicolson.
const double TAU = 0.005;                         // Time step.
const double T_FINAL = 2*TAU + 1e-4;              // Time interval length.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQ time step the mesh is unrefined.
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
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 5;                          // Maximum polynomial order allowed in hp-adaptivity
                                                  // had to be limited due to complicated integrals
const double ERR_STOP = 5.0;                      // Stopping criterion for hp-adaptivity
                                                  // (relative error between reference and coarse solution in percent)
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Newton's method.
const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 50;                   // Maximum allowed number of Newton iterations.

// Problem parameters.
const double H = 1;                               // Planck constant 6.626068e-34.
const double M = 1;                               // Mass of boson.
const double G = 1;                               // Coupling constant.
const double OMEGA = 1;                           // Frequency.


// Initial condition.
scalar init_cond(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-20*(x*x + y*y));
  dx = val * (-40.0*x);
  dy = val * (-40.0*y);
  return val;
}

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

// Weak forms.
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;

  // Convert initial condition into a Solution.
  Solution sln_prev_time;
  sln_prev_time.set_exact(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(J_euler), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_euler), HERMES_ANY, &sln_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_cranic), HERMES_ANY, &sln_prev_time);
  }

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec_coarse = new scalar[ndof];
  OGProjection::project_global(&space, &sln_prev_time, coeff_vec_coarse, matrix_solver);

  Solution init_proj;
  Solution::vector_to_solution(coeff_vec_coarse, &space, &init_proj);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  bool verbose = true;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

  // Cleanup after the Newton loop on the coarse mesh.
  delete matrix_coarse;
  delete rhs_coarse;
  delete solver_coarse;
  delete [] coeff_vec_coarse;
  
  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

      // Project on globally derefined mesh.
      info("Projecting previous fine mesh solution on derefined mesh.");
      OGProjection::project_global(&space, &sln_prev_time, &sln);
    }

    // Adaptivity loop:
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Space* ref_space = construct_refined_space(&space);

      scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
      DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
      }

      // Now we can deallocate the previous fine mesh.
      if(as > 1) delete ref_sln.get_mesh();

      // Newton's loop on the fine mesh.
      info("Solving on fine mesh:");
      if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
	  	        NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Store the result in ref_sln.
      Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
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
  }

  AbsFilter mag2(&sln);
  int success = 1;
  double eps = 1e-5;
  double val = std::abs(mag2.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.1, 0.1) xvel value = %lf", val);
  if (fabs(val - (0.65773)) > eps) {
    printf("Coordinate ( 0.1, 0.1) xvel value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.1, -0.1));
  info("Coordinate ( 0.1, -0.1) xvel value = %lf", val);
  if (fabs(val - (0.65773)) > eps) {
    printf("Coordinate ( 0.1, -0.1) xvel value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.2, 0.1));
  info("Coordinate ( 0.2, 0.1) xvel value = %lf", val);
  if (fabs(val - (0.39062)) > eps) {
    printf("Coordinate ( 0.2, 0.1) xvel value = %lf\n", val);
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
