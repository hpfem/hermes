#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;

// This test makes sure that example 40-trilinos-linear works correctly.

const int INIT_REF_NUM = 2;      // Number of initial uniform mesh refinements.
const int P_INIT = 3;            // Initial polynomial degree of all mesh elements.
const bool JFNK = false;         // true = Jacobian-free method (for NOX),
                                 // false = Newton (for NOX).
const bool PRECOND = true;       // Preconditioning by jacobian in case of JFNK (for NOX),
                                 // default ML preconditioner in case of Newton.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return x*x + y*y;
}

// Initial condition.
double init_cond(double x, double y, double &dx, double &dy)
{
	dx = 0;
	dy = 0;
	return 0;
}

// Exact solution.
double exact(double x, double y, double &dx, double &dy)
{
	dx = 2*x;
	dy = 2*y;
	return x*x +y*y;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **argv)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh* mesh = new Mesh();
  H2DReader mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
 
  // Create an H1 space with default shapeset.
  H1Space space(mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("---- Assembling by LinearProblem, solving by UMFpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation.
  WeakForm wf1;
  wf1.add_matrix_form(callback(bilinear_form));
  wf1.add_vector_form(callback(linear_form));

  // Initialize the linear problem.
  LinearProblem lp(&wf1, &space);

  // Select matrix solver.
  Matrix* mat; Vector* coeff_vec; CommonSolver* common_solver;
  init_matrix_solver(matrix_solver, ndof, mat, coeff_vec, common_solver);

  // Assemble stiffness matrix and rhs.
  lp.assemble(mat, coeff_vec);

  // Solve the matrix problem.
  if (!common_solver->solve(mat, coeff_vec)) error ("Matrix solver failed.\n");

  // Show the UMFpack solution.
  Solution sln_hermes(&space, coeff_vec);

  info("Coordinate (-0.6, -0.6) sln_hermes value = %lf", sln_hermes.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) sln_hermes value = %lf", sln_hermes.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) sln_hermes value = %lf", sln_hermes.get_pt_value( 0.0,  0.0));

  // CPU time needed by UMFpack.
  double umf_time = cpu_time.tick().last();

  // TRILINOS PART:
  info("---- Assembling by FeProblem, solving by NOX:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX to zero. Alternatively, you can obtain 
  // an initial vector by projecting init_cond() on the FE space, see below.
  coeff_vec->set_zero();

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();
  
  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  wf2.add_matrix_form(callback(jacobian_form), H2D_SYM);
  wf2.add_vector_form(callback(residual_form));

  // FIXME: The entire FeProblem should be removed
  // and functionality merged into LinearProblem.
  // Initialize FeProblem.
  FeProblem* fep = new FeProblem(&wf2, &space);

  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver* nox_solver = new NoxSolver(fep);
  nox_solver->set_init_sln(coeff_vec->get_c_array());

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver->set_precond(pc);
    else nox_solver->set_precond("ML");
  }

  // Assemble and solve using NOX.
  bool solved = nox_solver->solve();

  Solution sln_nox;
  if (solved)
  {
    double *coeffs = nox_solver->get_solution_vector();
    sln_nox.set_coeff_vector(&space, coeffs, ndof);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver->get_num_iters(), nox_solver->get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver->get_num_lin_iters(), nox_solver->get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Calculate errors.
  Solution ex;
  ex.set_exact(mesh, &exact);
  Adapt hp(&space, H2D_H1_NORM);
  hp.set_solutions(&sln_hermes, &ex);
  double err_est_rel_1 = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  info("Solution 1 (LinearProblem - UMFpack): exact H1 error: %g (time %g s)", err_est_rel_1, umf_time);
  hp.set_solutions(&sln_nox, &ex);
  double err_est_rel_2 = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  info("Solution 2 (FeProblem - NOX):  exact H1 error: %g (time %g + %g s)", 
    err_est_rel_2, proj_time, nox_time);

  info("Coordinate (-0.6, -0.6) sln_nox value = %lf", sln_nox.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) sln_nox value = %lf", sln_nox.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) sln_nox value = %lf", sln_nox.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) sln_nox value = %lf", sln_nox.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) sln_nox value = %lf", sln_nox.get_pt_value( 0.0,  0.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(sln_nox.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) sln_nox value is %g\n", sln_nox.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) sln_nox value is %g\n", sln_nox.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln_nox value is %g\n", sln_nox.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) sln_nox value is %g\n", sln_nox.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln_nox.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) sln_nox value is %g\n", sln_nox.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  if (fabs(sln_hermes.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) sln_hermes value is %g\n", sln_hermes.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) sln_hermes value is %g\n", sln_hermes.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln_hermes.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) sln_hermes value is %g\n", sln_hermes.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  delete nox_solver;
  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
