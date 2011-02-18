#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;

// This test makes sure that example 41-trilinos-nonlinear works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

const bool JFNK = false;                          // true = jacobian-free method,
                                                  // false = Newton.
const int PRECOND = 2;                            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                                  // in case of JFNK,
                                                  // Default ML proconditioner in case of Newton.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

// Initial condition.
double init_cond(double x, double y, double &dx, double &dy)
{
	dx = 0;
	dy = 0;
	return 0;
}

// Exact solution.
#include "exact_solution.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("Assembling by DiscreteProblem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation,
  WeakForm wf1;
  wf1.add_matrix_form(callback(jacobian_form_hermes), HERMES_NONSYM, HERMES_ANY);
  wf1.add_vector_form(callback(residual_form_hermes), HERMES_ANY);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp1(&wf1, &space, is_linear);

  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln_hermes;

  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)] ;
  Solution* sln_tmp = new Solution(&mesh, init_cond);
  OGProjection::project_global(&space, sln_tmp, coeff_vec, matrix_solver);
  delete sln_tmp;

  // Perform Newton's iteration.
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(&space);

    // Assemble the Jacobian matrix and residual vector.
    dp1.assemble(coeff_vec, matrix, rhs, false);

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
    if(!solver->solve())
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
    
    if (it >= NEWTON_MAX_ITER)
      error ("Newton method did not converge.");

    it++;
  }

  // Translate the resulting coefficient vector into the Solution sln_hermes.
  Solution::vector_to_solution(coeff_vec, &space, &sln_hermes);

  // Cleanup.
  delete(matrix);
  delete(rhs);
  delete(solver);

  // CPU time needed by UMFpack
  double umf_time = cpu_time.tick().last();

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);
 
  // TRILINOS PART:

  // Project the initial condition on the FE space.
  info("Projecting initial condition on the FE space.");
  // The NULL pointer means that we do not want the projection result as a Solution.
  sln_tmp = new Solution(&mesh, init_cond);
  OGProjection::project_global(&space, sln_tmp, coeff_vec, matrix_solver);
  delete sln_tmp;

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_matrix_form(callback(jacobian_form_nox), HERMES_SYM);
  if (JFNK && PRECOND == 2) wf2.add_matrix_form(callback(precond_form_nox), HERMES_SYM);
  wf2.add_vector_form(callback(residual_form_nox));

  // Initialize DiscreteProblem.
  DiscreteProblem dp2(&wf2, &space);

  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver nox_solver(&dp2);
  nox_solver.set_init_sln(coeff_vec);

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the matrix problem using NOX.
  info("Assembling by DiscreteProblem, solving by NOX.");
  bool solved = nox_solver.solve();
  Solution sln_nox;
  if (solved)
  {
    double *coeffs = nox_solver.get_solution();
    Solution::vector_to_solution(coeffs, &space, &sln_nox);
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
         nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
         nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else
    error("NOX failed.");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Calculate errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  Adapt adaptivity(&space);
  bool solutions_for_adapt = false;
  double err_est_rel_1 = adaptivity.calc_err_exact(&sln_hermes, &ex, solutions_for_adapt) * 100;
  info("Solution 1 (DiscreteProblem + %s): exact H1 error: %g (time %g [s])", MatrixSolverNames[matrix_solver].c_str(), err_est_rel_1, umf_time);
  double err_est_rel_2 = adaptivity.calc_err_exact(&sln_nox, &ex, solutions_for_adapt) * 100;
  info("Solution 2 (DiscreteProblem + NOX): exact H1 error: %g (time %g + %g = %g [s])", err_est_rel_2, proj_time, nox_time, proj_time+nox_time);

  info("Coordinate ( 0.6,  0.6) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.6,  0.6));
  info("Coordinate ( 0.4,  0.6) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.4,  0.6));
  info("Coordinate ( 0.4,  0.4) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.4,  0.4));
  info("Coordinate ( 0.6,  0.0) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.6,  0.0));
  info("Coordinate ( 0.5,  0.5) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.5,  0.5));

  info("Coordinate ( 0.6,  0.6) sln_nox value = %lf", sln_nox.get_pt_value( 0.6,  0.6));
  info("Coordinate ( 0.4,  0.6) sln_nox value = %lf", sln_nox.get_pt_value( 0.4,  0.6));
  info("Coordinate ( 0.4,  0.4) sln_nox value = %lf", sln_nox.get_pt_value( 0.4,  0.4));
  info("Coordinate ( 0.6,  0.0) sln_nox value = %lf", sln_nox.get_pt_value( 0.6,  0.0));
  info("Coordinate ( 0.5,  0.5) sln_nox value = %lf", sln_nox.get_pt_value( 0.5,  0.5));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(sln_nox.get_pt_value(0.6, 0.6) - 0.057600) > eps) {
    printf("Coordinate (0.6, 0.6) sln_nox value is %g\n", sln_nox.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.4, 0.6) - 0.057600) > eps) {
    printf("Coordinate ( 0.4, 0.6) sln_nox value is %g\n", sln_nox.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.4,  0.4) - 0.057600) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln_nox value is %g\n", sln_nox.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value(0.6,  0.0) - 0.000000) > eps) {
    printf("Coordinate (0.6,  0.0) sln_nox value is %g\n", sln_nox.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.5,  0.5) - 0.062500) > eps) {
    printf("Coordinate ( 0.5,  0.5) sln_nox value is %g\n", sln_nox.get_pt_value( 0.5,  0.5));
    success = 0;
  }

  if (fabs(sln_hermes.get_pt_value(0.6, 0.6) - 0.057600) > eps) {
    printf("Coordinate (0.6, 0.6) sln_hermes value is %g\n", sln_hermes.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.4, 0.6) - 0.057600) > eps) {
    printf("Coordinate ( 0.4, 0.6) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.4,  0.4) - 0.057600) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value(0.6,  0.0) - 0.000000) > eps) {
    printf("Coordinate (0.6,  0.0) sln_hermes value is %g\n", sln_hermes.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.5,  0.5) - 0.062500) > eps) {
    printf("Coordinate ( 0.5,  0.5) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.5,  0.5));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
