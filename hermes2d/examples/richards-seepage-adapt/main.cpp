#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  the time-dependent Richard's equation. The time discretization 
//  is backward Euler or Crank-Nicolson, and the Newton's method 
//  is applied to solve the nonlinear problem in each time step. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: rectangle (0, 8) x (0, 6.5).
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
// #define CONSTITUTIVE_GENUCHTEN

const int P_INIT = 1;                      // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 0;            // Number of initial mesh refinements towards the top edge.
const int TIME_INTEGRATION = 2;            // 1... implicit Euler, 2... Crank-Nicolson.

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
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Newton's method
const double NEWTON_TOL_COARSE = 0.0001;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.0005;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 50;            // Maximum allowed number of Newton iterations.

// Problem parameters.
const double TAU = 5e-3;                   // Time step.
const double STARTUP_TIME = 1.1e-2;        // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 5.0;                // Time interval length.
double TIME = 0;                           // Global time variable initialized with first time step.
double H_INIT = -9.5;                      // Initial pressure head.
double H_ELEVATION = 5.2;

double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
double K_S_4 = 1.061;

double ALPHA_1 = 0.01;
double ALPHA_3 = 0.005;
double ALPHA_2 = 0.01;
double ALPHA_4 = 0.05;

double THETA_R_1 = 0.1020;
double THETA_R_2 = 0.09849;
double THETA_R_3 = 0.08590;
double THETA_R_4 = 0.08590;

double THETA_S_1 = 0.4570;
double THETA_S_2 = 0.4510;
double THETA_S_3 = 0.4650;
double THETA_S_4 = 0.5650;

double N_1 = 1.982;
double N_2 = 1.632; 
double N_3 = 5.0;
double N_4 = 5.0;

double M_1 = 0.49546;
double M_2 = 0.38726;
double M_3 = 0.8;
double M_4 = 0.8;

double Q_MAX_VALUE = 0.07;         // Maximum value, used in function q_function(); 
double q_function() {
  if (STARTUP_TIME > TIME) return Q_MAX_VALUE * TIME / STARTUP_TIME;
  else return Q_MAX_VALUE;
}

double STORATIVITY = 0.05;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Material properties.
bool is_in_mat_1(double x, double y) {
  if (y >= -0.5) return true;
  else return false; 
}

bool is_in_mat_2(double x, double y) {
  if (y >= -1.0 && y < -0.5) return true;
  else return false; 
}

bool is_in_mat_4(double x, double y) {
  if (x >= 1.0 && x <= 3.0 && y >= -2.5 && y < -1.5) return true;
  else return false; 
}

bool is_in_mat_3(double x, double y) {
  if (!is_in_mat_1(x, y) && !is_in_mat_2(x, y) && !is_in_mat_4(x, y)) return true;
  else return false; 
}

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
int BDY_1 = 1;
int BDY_3 = 3;
int BDY_4 = 4;
int BDY_6 = 6;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == BDY_3) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Initial condition.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = 0;
  dy = -1;
  return -y + H_INIT;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  if (STARTUP_TIME > TIME) return -y + H_INIT + TIME/STARTUP_TIME*H_ELEVATION;
  else return -y + H_INIT + H_ELEVATION;
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
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(3, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution sln, ref_sln, sln_prev_time;
  
  // Adapt mesh to represent initial condition with given accuracy.
  info("Mesh adaptivity to an exact function:");
  // Initialize views.
  char title_init[200];
  sprintf(title_init, "Projection of initial condition");
  ScalarView* view_init = new ScalarView(title_init, 0, 0, 410, 300);
  sprintf(title_init, "Initial mesh");
  OrderView* ordview_init = new OrderView(title_init, 420, 0, 350, 300);
  view_init->fix_scale_width(80);
  int as = 1; bool done = false;
  do
  {
    // Setup space for the reference solution.
    Space *rspace = construct_refined_space(&space);

    // Assign the function f() to the fine mesh.
    ref_sln.set_exact(rspace->get_mesh(), init_cond);

    // Project the function f() on the coarse mesh.
    OGProjection::project_global(&space, &ref_sln, &sln_prev_time, matrix_solver);

    // Calculate element errors and total error estimate.
    Adapt adaptivity(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity.calc_err_est(&sln_prev_time, &ref_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    info("Step %d, ndof %d, proj_error %g%%", as, Space::get_num_dofs(&space), err_est_rel);

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      double to_be_processed = 0;
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, to_be_processed);

      if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

      view_init->show(&sln_prev_time);
      char title_init[100];
      sprintf(title_init, "Initial mesh, step %d", as);
      ordview_init->set_title(title_init);
      ordview_init->show(&space);
    }
    as++;
  }
  while (done == false);
  
  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain coefficient vector for Newton on coarse mesh.");
  scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(&space)];
  OGProjection::project_global(&space, init_cond, coeff_vec_coarse, matrix_solver);
  OGProjection::project_global(&space, &sln_prev_time, &sln, matrix_solver);

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_form_vol_euler, jac_form_vol_ord, HERMES_UNSYM, HERMES_ANY, 
                       &sln_prev_time);
    wf.add_matrix_form_surf(jac_form_surf_1_euler, jac_form_surf_1_ord, BDY_1);
    wf.add_matrix_form_surf(jac_form_surf_4_euler, jac_form_surf_4_ord, BDY_4);
    wf.add_matrix_form_surf(jac_form_surf_6_euler, jac_form_surf_6_ord, BDY_6);
    wf.add_vector_form(res_form_vol_euler, res_form_vol_ord, HERMES_ANY, 
                       &sln_prev_time);
    wf.add_vector_form_surf(res_form_surf_1_euler, res_form_surf_1_ord, BDY_1); 
    wf.add_vector_form_surf(res_form_surf_4_euler, res_form_surf_4_ord, BDY_4);
    wf.add_vector_form_surf(res_form_surf_6_euler, res_form_surf_6_ord, BDY_6);
  }
  else {
    wf.add_matrix_form(jac_form_vol_cranic, jac_form_vol_ord, HERMES_UNSYM, HERMES_ANY, 
                       &sln_prev_time);
    wf.add_matrix_form_surf(jac_form_surf_1_cranic, jac_form_surf_1_ord, BDY_1);
    wf.add_matrix_form_surf(jac_form_surf_4_cranic, jac_form_surf_4_ord, BDY_4);
    wf.add_matrix_form_surf(jac_form_surf_6_cranic, jac_form_surf_6_ord, BDY_6); 
    wf.add_vector_form(res_form_vol_cranic, res_form_vol_ord, HERMES_ANY, 
                       &sln_prev_time);
    wf.add_vector_form_surf(res_form_surf_1_cranic, res_form_surf_1_ord, BDY_1, 
			    &sln_prev_time);
    wf.add_vector_form_surf(res_form_surf_4_cranic, res_form_surf_4_ord, BDY_4, 
			    &sln_prev_time);
    wf.add_vector_form_surf(res_form_surf_6_cranic, res_form_surf_6_ord, BDY_6, 
			    &sln_prev_time);
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, graph_time_dof, graph_time_cpu;
 
  // Visualize the projection and mesh.
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(450, 0, 400, 350));
  view.show(&sln_prev_time);
  ordview.show(&space);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Time measurement.
    cpu_time.tick();

    // Updating current time.
    TIME = ts*TAU;
    info("---- Time step %d:", ts);

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
      FeProblem fep(&wf, ref_space, is_linear);

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
        fep.assemble(coeff_vec, matrix, rhs, false);

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
      info("Calculating error estimate."); 
      Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = true;
      
      // Calculate error estimate wrt. fine mesh solution.
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, space_err_est_rel: %g%%", 
        Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

      // Add entries to convergence graphs.
      graph_time_err_est.add_values(ts*TAU, err_est_rel);
      graph_time_err_est.save("time_error_est.dat");
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

    // Copy new time level solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
