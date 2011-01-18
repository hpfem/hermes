#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  a simple version of the time-dependent Richard's equation.
//  The time discretization is backward Euler, and the Newton's
//  method is applied to solve the nonlinear problem in each time 
//  step. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Known exact solution, see the file exact_solution.cpp.
//
//  Domain: square (0, 100)^2.
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
// Note: Exact solution makes sense for Gardner's relations only.
//#define CONSTITUTIVE_GENUCHTEN

const double TIME_INIT = 1e-3;                    // Initial time.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 8;                   // Number of initial mesh refinements towards the top edge.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double TAU = 0.001;                         // Time step.
const double T_FINAL = 1.0;                       // Time interval length.
const int TIME_INTEGRATION = 1;                   // 1... implicit Euler, 2... Crank-Nicolson.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Adaptivity
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
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
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const double ERR_STOP_INIT = 3.1;                 // Stopping criterion for the initial mesh adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;                   // Maximum allowed number of Newton iterations.

// Problem parameters.
double K_S = 20.464;
double ALPHA = 1e-3;
double THETA_R = 0;
double THETA_S = 0.45;
double H_R = -1000;
double A = 100;
double L = 100;
double STORATIVITY = 0.0;
double M = 0.6;
double N = 2.5;

// Current time.
double TIME = TIME_INIT;

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Exact solution (based on Fourier series)
#include "exact_solution.cpp"

// Boundary markers.
const int BDY_TOP = 2;
const int BDY_REST = 1;

// Boundary and initial conditions.
int Y_POWER = 5;
double bdy_cond(double x, double y, double& dx, double& dy) 
{
  return exact_sol(x, y, dx, dy);
  //dx = (100 - 2*x)/2.5 * pow(y/100, Y_POWER);
  //dy = x*(100 - x)/2.5 * pow(y/100, Y_POWER - 1) * 1./100;
  //return x*(100 - x)/2.5 * pow(y/100, Y_POWER) - 1000;
}
double init_cond(double x, double y, double& dx, double& dy) 
{
  return exact_sol(x, y, dx, dy);

  //dx = 0;
  //dy = 0;
  //return -1000;
}

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(double x, double y)
{
  double dx, dy;
  return bdy_cond(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_TOP, INIT_REF_NUM_BDY);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_TOP, BDY_REST));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(Hermes::vector<int>(BDY_TOP, BDY_REST), essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Create an H1 space for the initial coarse mesh solution.
  H1Space init_space(&basemesh, &bc_types, &bc_values, P_INIT);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and adaptivity.
  Solution sln_prev_time, sln, ref_sln;

  // Initialize views.
  char title_init[200];
  sprintf(title_init, "Projection of initial condition");
  ScalarView* view_init = new ScalarView(title_init, new WinGeom(0, 0, 410, 300));
  sprintf(title_init, "Initial mesh");
  OrderView* ordview_init = new OrderView(title_init, new WinGeom(420, 0, 350, 300));
  view_init->fix_scale_width(80);

  // Initialize sln_prev_time.
  // Note: only if adaptivity to initial condition is not done.
  sln_prev_time.set_exact(&basemesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, HERMES_NONSYM, HERMES_ANY, &sln_prev_time);
    wf.add_vector_form(res_euler, res_ord, HERMES_ANY, &sln_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, HERMES_NONSYM, HERMES_ANY, &sln_prev_time);
    wf.add_vector_form(res_cranic, res_ord, HERMES_ANY, &sln_prev_time);
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, graph_time_dof, graph_time_cpu;

  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain coefficient vector for Newton on coarse mesh.");
  scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(&space)];
  OGProjection::project_global(&space, init_cond, coeff_vec_coarse, matrix_solver);

  ScalarView view("Projection of initial condition", new WinGeom(0, 0, 410, 300));
  OrderView ordview("Initial mesh", new WinGeom(420, 0, 350, 300));
  view.fix_scale_width(80);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, &space, is_linear);

  // Set up the solver_coarse, matrix_coarse, and rhs_coarse according to the solver_coarse selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Perform Newton's iteration.
  info("Solving on coarse mesh.");
  bool verbose = true;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the actual solution. 
  Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

  // Clean up.
  delete [] coeff_vec_coarse;
  delete rhs_coarse;
  delete matrix_coarse;
  delete solver_coarse;

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Time measurement.
    cpu_time.tick();

    // Updating current time.
    TIME = ts*TAU;

    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      if (UNREF_LEVEL == 1) mesh.unrefine_all_elements();
      else mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }

    // Spatial adaptivity loop. Note; sln_prev_time must not be changed during 
    // spatial adaptivity.
    bool done = false;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Space* ref_space = construct_refined_space(&space);

      scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
     
      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1 && ts == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
        delete ref_sln.get_mesh();
      }

      // Initialize the FE problem.
      bool is_linear = false;
      DiscreteProblem dp(&wf, ref_space, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Perform Newton's iteration.
      info("Solving on fine mesh.");
      if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
          NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Translate the resulting coefficient vector into the actual solutions. 
      Solution::vector_to_solutions(coeff_vec, ref_space, &ref_sln);

      // Project the fine mesh solution on the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error calculation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt* adaptivity = new Adapt(&space);
      
      // Calculate error estimate wrt. fine mesh solution.
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Calculate error wrt. exact solution.
      ExactSolution exact(&mesh, exact_sol);
      bool solutions_for_adapt = false;
      double err_exact_rel = adaptivity->calc_err_exact(&sln, &exact, solutions_for_adapt) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d", 
	    Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
      info("space_err_est_rel: %g%%, space_err_exact_rel: %g%%", 
	    err_est_rel, err_exact_rel);

      // Add entries to convergence graphs.
      graph_time_err_est.add_values(ts*TAU, err_est_rel);
      graph_time_err_est.save("time_error_est.dat");
      graph_time_err_exact.add_values(ts*TAU, err_exact_rel);
      graph_time_err_exact.save("time_error_exact.dat");
      graph_time_dof.add_values(ts*TAU, Space::get_num_dofs(&space));
      graph_time_dof.save("time_dof.dat");
      graph_time_cpu.add_values(ts*TAU, cpu_time.accumulated());
      graph_time_cpu.save("time_cpu.dat");

      // If space_err_est too large, adapt the mesh.
      if (err_est_rel < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(&space) >= NDOF_STOP) {
          done = true;
          break;
        }

        as++;
      }

      // Cleanup.
      delete [] coeff_vec;
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      delete ref_space;
    }
    while (!done);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time level %d", ts);
    view.set_title(title);
    view.show(&sln);
    sprintf(title, "Mesh, time level %d", ts);
    ordview.set_title(title);
    ordview.show(&space);

    // Copy new time level solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
