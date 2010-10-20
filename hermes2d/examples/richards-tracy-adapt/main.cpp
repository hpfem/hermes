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

const double TIME_INIT = 1e-3;                    // Initial time.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 8;                   // Number of initial mesh refinements towards the top edge.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double TAU = 0.001;                         // Time step.
const double T_FINAL = 1.0;                       // Time interval length.
const int TIME_INTEGRATION = 1;                   // 1... implicit Euler, 2... Crank-Nicolson.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

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
                                           // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const double ERR_STOP_INIT = 3.1;          // Stopping criterion for the initial mesh adaptivity (rel. error tolerance between the
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
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and adaptivity.
  Solution u_prev_time, sln, ref_sln;

  // Initialize views.
  char title_init[200];
  sprintf(title_init, "Projection of initial condition");
  ScalarView* view_init = new ScalarView(title_init, 0, 0, 410, 300);
  sprintf(title_init, "Initial mesh");
  OrderView* ordview_init = new OrderView(title_init, 420, 0, 350, 300);
  view_init->fix_scale_width(80);

  /*
  // Adapt mesh to represent initial condition with given accuracy.
  info("Mesh adaptivity to an exact function:");
  int as = 1; bool done = false;
  do
  {
    // Setup space for the reference solution.
    Space *rspace = construct_refined_space(&space);

    // Assign the function f() to the fine mesh.
    ref_sln.set_exact(rspace->get_mesh(), init_cond);

    // Project the function f() on the coarse mesh.
    OGProjection::project_global(&space, &ref_sln, &u_prev_time, matrix_solver);

    // Calculate element errors and total error estimate.
    Adapt adaptivity(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity.calc_err_est(&u_prev_time, &ref_sln, solutions_for_adapt, 
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    info("Step %d, ndof %d, proj_error %g%%", as, Space::get_num_dofs(&space), err_est_rel);

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP_INIT) done = true;
    else {
      double to_be_processed = 0;
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, to_be_processed);

      if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

      view_init->show(&u_prev_time);
      char title_init[100];
      sprintf(title_init, "Initial mesh, step %d", as);
      ordview_init->set_title(title_init);
      ordview_init->show(&space);
    }
    as++;
  }
  while (done == false);
  */

  // Initialize u_prev_time.
  // Note: only if adaptivity to initial condition is not done.
  u_prev_time.set_exact(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, HERMES_UNSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_euler, res_ord, HERMES_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, HERMES_UNSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_cranic, res_ord, HERMES_ANY, &u_prev_time);
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, graph_time_dof, graph_time_cpu;

  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain coefficient vector for Newton on coarse mesh.");
  scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(&space)];
  OGProjection::project_global(&space, init_cond, coeff_vec_coarse, matrix_solver);

  ScalarView view("Projection of initial condition", 0, 0, 410, 300);
  OrderView ordview("Initial mesh", 420, 0, 350, 300);
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
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(&space);

    // Assemble the Jacobian matrix_coarse and residual vector.
    dp_coarse.assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse, false);

    // Multiply the residual vector with -1 since the matrix_coarse 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs_coarse->set(i, -rhs_coarse->get(i));
    
    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs_coarse);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(&space), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if (res_l2_norm < NEWTON_TOL_COARSE || it > NEWTON_MAX_ITER) break;

    // Solve the linear system.
    if(!solver_coarse->solve())
      error ("matrix_coarse solver_coarse failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];
    
    if (it >= NEWTON_MAX_ITER)
      error ("Newton method did not converge.");

    it++;
  }

  // Translate the resulting coefficient vector into the actual solutions. 
  Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);
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

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

      // Project fine mesh solution on the globally derefined mesh.
      info("Projecting fine mesh solution on globally derefined mesh.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);
    }

    // Adaptivity loop (in space):
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
      }

      // Initialize the FE problem.
      bool is_linear = false;
      DiscreteProblem dp(&wf, ref_space, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Perform Newton's iteration.
      int it = 1;
      while (1)
      {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(ref_space);

        // Assemble the Jacobian matrix and residual vector.
        dp.assemble(coeff_vec, matrix, rhs, false);

        // Multiply the residual vector with -1 since the matrix 
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
        
        // Calculate the l2-norm of residual vector.
        double res_l2_norm = get_l2_norm(rhs);

        // Info for user.
        info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(ref_space), res_l2_norm);

        // If l2 norm of the residual vector is within tolerance, or the maximum number 
        // of iteration has been reached, then quit.
        if (res_l2_norm < NEWTON_TOL_FINE || it > NEWTON_MAX_ITER) break;

        // Solve the linear system.
        if(!solver->solve())
          error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
        
        if (it >= NEWTON_MAX_ITER)
          error ("Newton method did not converge.");

        it++;
      }
      
      // Translate the resulting coefficient vector into the actual solutions. 
      Solution::vector_to_solutions(coeff_vec, ref_space, &ref_sln);

      // Project the fine mesh solution on the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error calculation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = true;
      
      // Calculate error estimate wrt. fine mesh solution.
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                           HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

      // Calculate error wrt. exact solution.
      ExactSolution exact(&mesh, exact_sol);
      solutions_for_adapt = false;
      double err_exact_rel = adaptivity->calc_err_exact(&sln, &exact, solutions_for_adapt, 
                             HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

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
      delete ref_space->get_mesh();
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

    // Copy new time level solution into u_prev_time.
    u_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
