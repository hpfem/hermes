#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;

// This test makes sure that example 41-trilinos-nonlinear works correctly.

const int INIT_REF_NUM = 2;       // Number of initial uniform mesh refinements.
const int P_INIT = 3;             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations.

const bool JFNK = false;          // true = jacobian-free method,
                                  // false = Newton.
const int PRECOND = 2;            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                  // in case of JFNK,
                                  // Default ML proconditioner in case of Newton.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

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

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, NULL, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("Assembling by DiscreteProblem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation,
  WeakForm wf1;
  wf1.add_matrix_form(callback(jacobian_form_hermes), H2D_UNSYM, H2D_ANY);
  wf1.add_vector_form(callback(residual_form_hermes), H2D_ANY);

  // Initialize NonlinSystem,
  DiscreteProblem dp(&wf1, &space);

  // Select matrix solver.
  Matrix* mat; Vector* coeff_vec; CommonSolver* common_solver;
  init_matrix_solver(matrix_solver, ndof, mat, coeff_vec, common_solver);

  // Project the initial condition on the FE space.
  info("Projecting initial condition on the FE space.");
  // The NULL pointer means that we do not want the projection result as a Solution.
  Solution* sln_tmp = new Solution(&mesh, init_cond);
  project_global(&space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
  delete sln_tmp;

  // Perform Newton's iteration,
  info("Performing Newton's method.");
  bool verbose = true;
  if (!solve_newton(&space, &wf1, coeff_vec, matrix_solver, 
		    NEWTON_TOL, NEWTON_MAX_ITER, verbose)) {
    error("Newton's method did not converge.");
  };

  // Store the solution in "sln_hermes".
  Solution sln_hermes(&space, coeff_vec);

  // CPU time needed by UMFpack
  double umf_time = cpu_time.tick().last();

  info("Assembling by FeProblem, solving by NOX:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);
 
  // TRILINOS PART:

  // Project the initial condition on the FE space.
  info("Projecting initial condition on the FE space.");
  // The NULL pointer means that we do not want the projection result as a Solution.
  sln_tmp = new Solution(&mesh, init_cond);
  project_global(&space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
  delete sln_tmp;

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_matrix_form(callback(jacobian_form_nox), H2D_SYM);
  if (JFNK && PRECOND == 2) wf2.add_matrix_form(callback(precond_form_nox), H2D_SYM);
  wf2.add_vector_form(callback(residual_form_nox));

  // Initialize FeProblem.
  FeProblem fep(&wf2, &space);

  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver nox_solver(&fep);
  nox_solver.set_init_sln(coeff_vec->get_c_array());

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the matrix problem using NOX.
  info("Assembling by FeProblem, solving by NOX.");
  bool solved = nox_solver.solve();
  Solution sln_nox;
  if (solved)
  {
    double *s = nox_solver.get_solution_vector();
    AVector *tmp_vector = new AVector(ndof);
    tmp_vector->set_c_array(s, ndof);
    sln_nox.set_fe_solution(&space, tmp_vector);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
         nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
         nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else
    error("NOX failed.");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Calculate exact errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  Adapt hp(&space, H2D_H1_NORM);
  hp.set_solutions(&sln_hermes, &ex);
  double err_est_rel_1 = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  hp.set_solutions(&sln_nox, &ex);
  double err_est_rel_2 = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  info("Solution 1 (DiscreteProblem + UMFpack): exact H1 error: %g (time %g s)", err_est_rel_1, umf_time);
  info("Solution 2 (FeProblem + NOX):  exact H1 error: %g (time %g + %g s)", err_est_rel_2, proj_time, nox_time);

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
