#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
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

const double TIME_INIT = 1e-3;             // Initial time.
const int INIT_REF_NUM = 0;                // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 8;            // Number of initial mesh refinements towards the top edge.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const double TAU = 0.001;                  // Time step.
const double T_FINAL = 1.0;                // Time interval length.
const int TIME_INTEGRATION = 1;            // 1... implicit Euler, 2... Crank-Nicolson.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;    // Predefined list of element refinement candidates. Possible values are
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
                                           // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;              // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

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

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

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
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
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
  mesh.refine_towards_boundary(2, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution u_prev_time;

  // Adapt mesh to represent initial condition with given accuracy.
  bool verbose = true;        // Report results.
  double err_stop_temp = 5.0;  // 
  adapt_to_exact_function(space, H2D_H1_NORM, init_cond, &selector, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY, err_stop_temp, NDOF_STOP, 
                          verbose, &u_prev_time);

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, H2D_UNSYM, H2D_ANY, &u_prev_time);
    wf.add_vector_form(res_euler, res_ord, H2D_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, H2D_UNSYM, H2D_ANY, &u_prev_time);
    wf.add_vector_form(res_cranic, res_ord, H2D_ANY, &u_prev_time);
  }

  // Initialize matrix solver.
  Matrix* mat; Vector* coeff_vec; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, coeff_vec, solver);

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, graph_time_dof, graph_time_cpu;

  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  Solution* sln_tmp = new Solution(&mesh, init_cond);
  // The NULL means that we do not want the result as a Solution.
  project_global(space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
  delete sln_tmp;

  // Initialize views.
  ScalarView view("Projection of initial condition", 0, 0, 410, 300);
  OrderView ordview("Initial mesh", 420, 0, 350, 300);
  view.fix_scale_width(80);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");
  if (!solve_newton(space, &wf, coeff_vec, matrix_solver, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the result in sln.
  Solution sln, ref_sln;
  sln.set_fe_solution(space, coeff_vec);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Time measurement.
    cpu_time.tick();

    // Updating current time.
    TIME = ts*TAU;

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space->set_uniform_order(P_INIT);

      // Project fine mesh solution on the globally derefined mesh.
      info("---- Time step %d:", ts);
      info("Projecting fine mesh solution on globally derefined mesh.");
      // The NULL means that we do not want the coefficient vector.
      project_global(space, H2D_H1_NORM, &ref_sln, &sln, NULL);
    }

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est_rel, space_err_exact_rel;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Mesh *ref_mesh = new Mesh();
      ref_mesh->copy(space->get_mesh());
      ref_mesh->refine_all_elements();
      Space* ref_space = space->dup(ref_mesh);
      int order_increase = 1;
      ref_space->copy_orders(space, order_increase);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1 && ts == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        // The NULL means that we do not want the result as a Solution.
        project_global(ref_space, H2D_H1_NORM, &sln, NULL, coeff_vec);
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        // The NULL means that we do not want the result as a Solution.
        project_global(ref_space, H2D_H1_NORM, &ref_sln, NULL, coeff_vec);
      }

      // Newton's method on fine mesh
      info("Solving on fine mesh.");
      if (!solve_newton(ref_space, &wf, coeff_vec, matrix_solver, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in ref_sln.
      ref_sln.set_fe_solution(ref_space, coeff_vec);

      // Calculate error estimate wrt. fine mesh solution.
      info("Calculating error (est).");
      Adapt hp(space, H2D_H1_NORM);
      hp.set_solutions(&sln, &ref_sln);
      double space_err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

      // Calculate error wrt. exact solution.
      info("Calculating error (exact).");
      ExactSolution exact(&mesh, exact_sol);
      double err_exact_abs = calc_abs_error(&sln, &exact, H2D_H1_NORM);
      double norm_exact = calc_norm(&exact, H2D_H1_NORM);
      space_err_exact_rel = err_exact_abs / norm_exact * 100;    

      info("ndof_coarse: %d, ndof_fine: %d, space_err_est_rel: %g%%, space_err_exact_rel: %g%%", 
	   get_num_dofs(space), get_num_dofs(ref_space), space_err_est_rel, space_err_exact_rel);

      // If space_err_est too large, adapt the mesh.
      if (space_err_est_rel < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (get_num_dofs(space) >= NDOF_STOP) {
          done = true;
          break;
        }

        // Project the fine mesh solution on the new coarse mesh.
        info("Projecting fine mesh solution on coarse mesh for error calculation.");
        // The NULL means that we do not want the coefficient vector.
        project_global(space, H2D_H1_NORM, &ref_sln, &sln, NULL);

        as++;
      }
    }
    while (!done);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time level %d", ts);
    view.set_title(title);
    view.show(&sln);
    sprintf(title, "Mesh, time level %d", ts);
    ordview.set_title(title);
    ordview.show(space);

    // Add entries to convergence graphs.
    graph_time_err_est.add_values(ts*TAU, space_err_est_rel);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_err_exact.add_values(ts*TAU, space_err_exact_rel);
    graph_time_err_exact.save("time_error_exact.dat");
    graph_time_dof.add_values(ts*TAU, get_num_dofs(space));
    graph_time_dof.save("time_dof.dat");
    graph_time_cpu.add_values(ts*TAU, cpu_time.accumulated());
    graph_time_cpu.save("time_cpu.dat");

    // Copy new time level solution into u_prev_time.
    u_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
