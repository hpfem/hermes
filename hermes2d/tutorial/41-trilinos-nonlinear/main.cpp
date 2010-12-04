#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
//  compares performance of the Newton's method in Hermes (assembling via the DiscreteProblem 
//  class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX 
//  solver (using the Hermes DiscreteProblem class to assemble discrete problems).
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

using namespace Teuchos;

const int INIT_REF_NUM = 2;       // Number of initial uniform mesh refinements.
const int P_INIT = 3;             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations.

const bool JFNK = false;          // true = jacobian-free method,
                                  // false = Newton.
const int PRECOND = 2;            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                  // in case of JFNK,
                                  // Default ML proconditioner in case of Newton.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_AZTECOO,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

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
  bc_types.add_bc_dirichlet(Hermes::Tuple<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, NULL, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("Assembling by DiscreteProblem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation,
  WeakForm wf1;
  wf1.add_matrix_form(callback(jacobian_form_hermes), HERMES_UNSYM, HERMES_ANY);
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
    for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
    
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
  Solution sln_nox;
  if (nox_solver.solve())
  {
    Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln_nox);
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
  Adapt adaptivity(&space, HERMES_H1_NORM);
  bool solutions_for_adapt = false;
  double err_est_rel_1 = adaptivity.calc_err_exact(&sln_hermes, &ex, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
  info("Solution 1 (DiscreteProblem + %s): exact H1 error: %g (time %g [s])", MatrixSolverNames[matrix_solver].c_str(), err_est_rel_1, umf_time);
  double err_est_rel_2 = adaptivity.calc_err_exact(&sln_nox, &ex, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
  info("Solution 2 (DiscreteProblem + NOX): exact H1 error: %g (time %g + %g = %g [s])", err_est_rel_2, proj_time, nox_time, proj_time+nox_time);

  // Show both solutions.
  ScalarView view1("Solution 1", new WinGeom(0, 0, 500, 400));
  view1.show(&sln_hermes);
  ScalarView view2("Solution 2", new WinGeom(510, 0, 500, 400));
  view2.show(&sln_nox);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
