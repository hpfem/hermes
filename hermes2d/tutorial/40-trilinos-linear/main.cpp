#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

//  The purpose of this example is to show how to use Trilinos
//  for linear PDE problems. It compares performance of the LinearProblem 
//  class in Hermes using the UMFpack matrix solver with the performance
//  of the Trilinos NOX solver (using Newton's method or JFNK, with or 
//  without preconditioning).
//
//  PDE: Poisson equation.
//
//  Domain: Square (-1, 1)^2.
//
//  BC: Nonhomogeneous Dirichlet, see the function essential_bc_values() below.
//
//  Exact solution: x*x + y*y.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 6;      // Number of initial uniform mesh refinements.
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
  ScalarView sv("Solution", 0, 0, 440, 350);
  sv.show(&sln_hermes);

  // CPU time needed by UMFpack.
  double umf_time = cpu_time.tick().last();

  // TRILINOS PART:

  info("---- Assembling by FeProblem, solving by NOX:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX to zero. Alternatively, you can obtain 
  // an initial vector by projecting init_cond() on the FE space, see below.
  coeff_vec->set_zero();

  /* // Generate an initial vector for NOX by projecting init_cond().
  // Project the initial condition on the FE space. 
  info("Projecting initial solution on the FE mesh.");
  // The NULL pointer means that we do not want the projection result as a Solution.
  Solution* sln_tmp = new Solution(mesh, init_cond);
  project_global(&space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
  delete sln_tmp;
  */

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
  MlPrecond pc("sa");
  if (PRECOND)
  {
    if (JFNK) nox_solver->set_precond(&pc);
    else nox_solver->set_precond("ML");
  }

  // Assemble and solve using NOX.
  bool solved = nox_solver->solve();

  Solution sln_nox;
  if (solved)
  {
    double *s = nox_solver->get_solution_vector();
    Vector *tmp_vector = new AVector(ndof);
    tmp_vector->set_c_array(s, ndof);
    sln_nox.set_fe_solution(&space, tmp_vector);
    delete tmp_vector;
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver->get_num_iters(), nox_solver->get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver->get_num_lin_iters(), nox_solver->get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Show the NOX solution.
  ScalarView view2("Solution 2", 450, 0, 440, 350);
  view2.set_min_max_range(0, 2);
  view2.show(&sln_nox);

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

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
