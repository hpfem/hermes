#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_sdirk-22 Benchmarks/sdirk-22
 *  \{
 *  \brief This test makes sure that the benchmark "sdirk-22" works correctly.
 *
 *  \section s_params Parameters
 *   - INIT_GLOB_REF_NUM = 0
 *   - INIT_BDY_REF_NUM = 0
 *   - P_INIT = 2
 *   - TAU = 0.1
 *   - T_FINAL = 1.0
 *   - NEWTON_TOL = 1e-10
 *   - NEWTON_MAX_ITER = 100
 *   - matrix_solver = SOLVER_UMFPACK
 *
 *  \section s_res Results
 *   - DOFs: 529
 *   - Adaptivity steps: 10 
 */


const int INIT_GLOB_REF_NUM = 0;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double TAU = 0.1;                                 // Time step.
const double T_FINAL = 1.0;                       // Time interval length.
const double NEWTON_TOL = 1e-10;                  // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double LX = 1.0;                            // Domain size in the x and y dimensions.
const double LY = 1.0;                            // Note: must be in agreement with mesh file.
const double OMEGA = 1.0;

double TIME = 0.0;

const double SIGMA = 1.0;
const double ALPHA = 0.0;
const double GAMMA = 1 - 1/sqrt((double)2);
enum TimeDiscretization {IE, SDIRK};
TimeDiscretization method = SDIRK;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1.0 + pow(u, ALPHA);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return ALPHA*pow(u, ALPHA-1);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = 0.0;
  dy = 0.0;
  return 0.0;
}

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Heat sources (forcing term in accordance with exact solution).
template<typename Real>
Real heat_src(Real x, Real y, double t)
{
  return -ALPHA*pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LY-y)-exp(OMEGA*t)*(LY+y)*(LX-x)*(LY-y),2.0)*pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)*(LY-y),ALPHA-1.0)-ALPHA*pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)-exp(OMEGA*t)*(LX+x)*(LX-x)*(LY-y),2.0)*pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)*(LY-y),ALPHA-1.0)+exp(OMEGA*t)*(pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)*(LY-y),ALPHA)+1.0)*(LX+x)*(LX-x)*2.0+exp(OMEGA*t)*(pow(exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)*(LY-y),ALPHA)+1.0)*(LY+y)*(LY-y)*2.0+OMEGA*exp(OMEGA*t)*(LX+x)*(LY+y)*(LX-x)*(LY-y);
}

// Exact solution.
scalar exact_solution(double x, double y, double& dx, double& dy)
{
  dx = (-2*x)*(LY+y)*(LY-y)*exp(OMEGA*TIME);
  dy = (LY+x)*(LX-x)*(-2*y)*exp(OMEGA*TIME);
  return (LX+x)*(LX-x)*(LY+y)*(LY-y)*exp(OMEGA*TIME);
}

// Initial condition.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return exact_solution(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[]) 
{
  // This is a hack I used to run the code a dozen of times when plotting convergence graphs.
  if (argc > 1) {
    if (argv[1][0] == 'e') method = IE;
    else if (argv[1][0] == 's') method = SDIRK;
    else error("what are you doing?");
  }
  if (argc > 2) {
    TAU = std::atof(argv[2]);
  }

  // This is important to make sure we compare solution at exact same point in time when studying convergence.
  int N_STEP = std::ceil(T_FINAL / TAU);
  if (fabs(T_FINAL - N_STEP * TAU) > 1e-10) {
    error("bad choice of TAU");
  }

  info("t_final = %g, tau = %g, n = %i", T_FINAL, TAU, N_STEP);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);
  
  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  Solution u_prev_time(&mesh, exact_solution);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec, matrix_solver);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  if (method == IE) {
    info("IMPLICIT EULER METHOD");
    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(jac), HERMES_NONSYM, HERMES_ANY);
    wf.add_vector_form(callback(res), HERMES_ANY, &u_prev_time);

    // Initialize the FE problem. 
    bool is_linear = false; 
    DiscreteProblem dp(&wf, &space, is_linear); 

    // Time stepping loop:
    int ts = 0;
    do {
      info("---- Time step %d, t = %g s.", ++ts, TIME);

      info("We are computing solution at next time step TIME+TAU = %g s.", TIME+TAU);

      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solution.
      Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

      // Update time.
      TIME = TIME + TAU;

      // Compute exact error.
      Solution exact_sln(&mesh, exact_solution);
      double exact_h1_error = calc_abs_error(&u_prev_time, &exact_sln, HERMES_H1_NORM);
      info("TIME: %g s.", TIME);
      info("Exact error in l2-norm: %g.", exact_h1_error);
    } 
    while (ts < N_STEP);

    // Cleanup.
    delete [] coeff_vec;
    delete matrix;
    delete rhs;
    delete solver;

    // Hack to extract error at final time for convergence graph.
    Solution citrouille(&mesh, exact_solution);
    info("IE: tau %g, abs_error %g.", TAU, calc_abs_error(&u_prev_time, &citrouille, HERMES_H1_NORM));
  }
  else if (method == SDIRK) {
    info("SDIRK22");

    Solution Y1(&mesh, init_cond);
    Solution Y2(&mesh, init_cond);

    scalar* coeff_vec1 = new scalar[ndof];
    OGProjection::project_global(&space, &u_prev_time, coeff_vec1, matrix_solver);
    scalar* coeff_vec2 = new scalar[ndof];
    OGProjection::project_global(&space, &u_prev_time, coeff_vec2, matrix_solver);

    WeakForm wf1;
    wf1.add_matrix_form(callback(jac_Y), HERMES_NONSYM, HERMES_ANY);
    wf1.add_vector_form(callback(res_Y1), HERMES_ANY, Hermes::vector<MeshFunction*>(&u_prev_time));
    WeakForm wf2;
    wf2.add_matrix_form(callback(jac_Y), HERMES_NONSYM, HERMES_ANY);
    wf2.add_vector_form(callback(res_Y2), HERMES_ANY, Hermes::vector<MeshFunction*>(&u_prev_time, &Y1));

    // Initialize the FE problem. 
    bool is_linear = false;
    DiscreteProblem dp1(&wf1, &space, is_linear);
    DiscreteProblem dp2(&wf2, &space, is_linear);

    double current_time = 0.0; int ts = 0;
    do {
      info("---- Time step %d, t = %g s.", ++ts, current_time);

      // Compute Y1.
      info("Compute Y1 at t = %g s.", TIME+GAMMA*TAU);

      // Perform Newton's iteration for Y1.
      int it = 1;
      while (1) {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(&space);

        // Assemble the Jacobian matrix and residual vector.
        dp1.assemble(coeff_vec1, matrix, rhs, false);

        // Multiply the residual vector with -1 since the matrix 
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        rhs->change_sign();


        // Calculate the l2-norm of residual vector.
        double res_l2_norm = get_l2_norm(rhs);

        // Info for user.
        info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(&space), res_l2_norm);

        // If l2 norm of the residual vector is within tolerance, or the maximum number 
        // of iteration has been reached, then quit.
        if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

        // Solve the linear system.
        if(!solver->solve()) error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++) coeff_vec1[i] += solver->get_solution()[i];

        if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");

        it++;
      }

      // Store Y1.
      Solution::vector_to_solution(coeff_vec1, &space, &Y1);

      // Compute Y2.
      info("Compute Y2 at t = %g s.", TIME+TAU);

      // Perform Newton's iteration.
      while (1) {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(&space);

        // Assemble the Jacobian matrix and residual vector.
        dp2.assemble(coeff_vec2, matrix, rhs, false);

        // Multiply the residual vector with -1 since the matrix 
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        rhs->change_sign();

        // Calculate the l2-norm of residual vector.
        double res_l2_norm = get_l2_norm(rhs);

        // Info for user.
        info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(&space), res_l2_norm);

        // If l2 norm of the residual vector is within tolerance, or the maximum number 
        // of iteration has been reached, then quit.
        if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

        // Solve the linear system.
        if(!solver->solve()) error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++) coeff_vec2[i] += solver->get_solution()[i];

        if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");

        it++;
      }

      // Store Y2.
      Solution::vector_to_solution(coeff_vec2, &space, &Y2);

      // Update previous time level solution.
      u_prev_time = Y2;

      // Update time.
      TIME = TIME + TAU;

      // Compute exact error.
      Solution exact_sln(&mesh, exact_solution);
      double exact_h1_error = calc_abs_error(&u_prev_time, &exact_sln, HERMES_H1_NORM);
      info("TIME: %g s.", TIME);
      info("Exact error in H1-norm: %g.", exact_h1_error);
    } while (ts < N_STEP);

    // Hack to extract error at final time for convergence graph.
    Solution citrouille(&mesh, exact_solution);
    info("SDIRK: tau %g, abs_error %g.", TAU, calc_abs_error(&u_prev_time, &citrouille, HERMES_H1_NORM));
  }

  info("Coordinate ( 0.0, 0.0) u_prev_time value = %lf", u_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 0.3, 0.0) u_prev_time value = %lf", u_prev_time.get_pt_value(0.3, 0.0));
  info("Coordinate ( 0.5, 0.0) u_prev_time value = %lf", u_prev_time.get_pt_value(0.5, 0.0));
  info("Coordinate ( 0.7, 0.0) u_prev_time value = %lf", u_prev_time.get_pt_value(0.7, 0.0));
  info("Coordinate ( 0.9, 0.0) u_prev_time value = %lf", u_prev_time.get_pt_value(0.9, 0.0));

  double values[5];
  if (method == IE) {
    values[0] = 2.730284;
    values[1] = 2.484559; 
    values[2] = 2.047713;
    values[3] = 1.392445;
    values[4] = 0.518754;
  }
  else {
    values[0] = 2.717647;
    values[1] = 2.473058; 
    values[2] = 2.038235;
    values[3] = 1.386000;
    values[4] = 0.516353;
  }  

  int success = 1;
  double eps = 1e-5;
  if (fabs(u_prev_time.get_pt_value(0.0, 0.0) - (values[0])) > eps) {
    printf("Coordinate ( 0.0, 0.0) u_prev_time value = %lf\n", u_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  if (fabs(u_prev_time.get_pt_value(0.3, 0.0) - (values[1])) > eps) {
    printf("Coordinate ( 0.3, 0.0) u_prev_time value = %lf\n", u_prev_time.get_pt_value(0.3, 0.0));
    success = 0;
  }
  if (fabs(u_prev_time.get_pt_value(0.5, 0.0) - (values[2])) > eps) {
    printf("Coordinate ( 0.5, 0.0) u_prev_time value = %lf\n", u_prev_time.get_pt_value(0.5, 0.0));
    success = 0;
  }
  if (fabs(u_prev_time.get_pt_value(0.7, 0.0) - (values[3])) > eps) {
    printf("Coordinate ( 0.7, 0.0) u_prev_time value = %lf\n", u_prev_time.get_pt_value(0.7, 0.0));
    success = 0;
  }
  if (fabs(u_prev_time.get_pt_value(0.9, 0.0) - (values[4])) > eps) {
    printf("Coordinate ( 0.9, 0.0) u_prev_time value = %lf\n", u_prev_time.get_pt_value(0.9, 0.0));
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
