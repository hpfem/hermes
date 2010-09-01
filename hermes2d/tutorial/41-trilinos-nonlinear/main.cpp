#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
//  compares performance of The Newton's method in the NonlinSystem class in Hermes (using the 
//  Umfpack solver) with the performance of the Trilinos/NOX solver (using the Newton's method 
//  or JFNK, and with or without preconditioning).
//
//  PDE:  - \nabla (k \nabla u) = f
//  k = (1 + sqr(u_x) + sqr(u_y))^{-0.5}
//
//  Domain: Unit square.
//
//  BC: zero Dirichlet.
//
//  Exact solution: (x - x*x) * (y - y*y).
//
//  Initial guess for the Newton's method: zero function.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 5;       // Number of initial uniform mesh refinements.
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
  MlPrecond pc("sa");
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(&pc);
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

  // Show both solutions.
  ScalarView view1("Solution 1", 0, 0, 500, 400);
  view1.show(&sln_hermes);
  ScalarView view2("Solution 2", 510, 0, 500, 400);
  view2.show(&sln_nox);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
