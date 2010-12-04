#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear time-dependent PDE discretized implicitly in time
//  using implicit Euler or Crank-Nicolson.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f.
//
//  Domain: square (-10,10)^2.
//
//  BC:  Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;                // Number of initial uniform mesh refinements.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const int TIME_DISCR = 1;                  // 1 for implicit Euler, 2 for Crank-Nicolson.
const double TAU = 0.5;                    // Time step. Note: The Crank-Nicolson method is known 
                                           // to have problems with large time steps on coarse meshes. 
                                           // Do not use it here with TAU > 0.1.
const double T_FINAL = 2.0;                // Time interval length.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
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
const double ERR_STOP = 1.0;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_AZTECOO,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;           // Maximum allowed number of Newton iterations.

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
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary markers.
const int BDY_ESSENTIAL = 1;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
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
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_ESSENTIAL);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;

  // Convert initial condition into a Solution.
  Solution sln_prev_time;
  sln_prev_time.set_exact(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(J_euler), HERMES_UNSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_euler), HERMES_ANY, &sln_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), HERMES_UNSYM, HERMES_ANY);
    wf.add_vector_form(callback(F_cranic), HERMES_ANY, &sln_prev_time);
  }

  // Initialize the discrete problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, &space, is_linear);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(445, 0, 410, 350));
  view.show(&sln_prev_time);
  ordview.show(&space);
  
  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
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
      if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
          NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");
      Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

      // Cleanup after the Newton loop on the coarse mesh.
      delete matrix_coarse;
      delete rhs_coarse;
      delete solver_coarse;
      delete [] coeff_vec_coarse;
    }

    // Adaptivity loop:
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = construct_refined_space(&space);

      // Initialize matrix solver.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
      scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];

      // Initialize discrete problem on reference mesh.
      DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);

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
      if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
	  	        NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Store the result in ref_sln.
      Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = true;
      double err_est_rel_total = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

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

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g", ts*TAU);
    view.set_title(title);
    view.show_mesh(false);
    view.show(&sln);
    sprintf(title, "Mesh, time %g", ts*TAU);
    ordview.set_title(title);
    ordview.show(&space);

    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
