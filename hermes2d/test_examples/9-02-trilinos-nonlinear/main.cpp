#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Teuchos;

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
//  compares performance of the Newton's method in Hermes (assembling via the DiscreteProblem 
//  class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX 
//  solver (using the Hermes DiscreteProblem class to assemble discrete problems).
//
//  PDE:  - \nabla (k \nabla u) - f = 0
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

const int INIT_REF_NUM = 5;                       // Number of initial uniform mesh refinements.
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
// NOX parameters.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  // NOX error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-8;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-8;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc("Bdy", 0.0);
  EssentialBCs bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("Assembling by DiscreteProblem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation,
  CustomWeakForm wf1;

  // Initialize the discrete problem.
  DiscreteProblem dp1(&wf1, &space);
  
  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln1;

  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  // We can start with a zero vector.
  memset(coeff_vec, 0, ndof * sizeof(double));
  // Or we can project the initial condition to obtain the initial
  // coefficient vector.
  //info("Projecting to obtain initial vector for the Newton's method.");
  //CustomInitialSolution sln_tmp(&mesh);
  //OGProjection::project_global(&space, &sln_tmp, coeff_vec, matrix_solver);

  // Perform Newton's iteration.
  bool verbose = true;
  bool jacobian_changed = true;
  if (!hermes2d.solve_newton(coeff_vec, &dp1, solver, matrix, rhs, jacobian_changed,
      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln1.
  Solution::vector_to_solution(coeff_vec, &space, &sln1);

  // Cleanup.
  delete(matrix);
  delete(rhs);
  delete(solver);

  // CPU time needed by UMFpack
  double time1 = cpu_time.tick().last();

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);
 
  // Show UMFPACK solution.
  ScalarView view1("Solution 1", new WinGeom(0, 0, 500, 400));
  view1.show(&sln1);

  // Calculate error.
  CustomExactSolution ex(&mesh);
  double rel_err_1 = hermes2d.calc_rel_error(&sln1, &ex, HERMES_H1_NORM) * 100;
  info("Solution 1 (%s):  exact H1 error: %g%% (time %g s)", MatrixSolverNames[matrix_solver].c_str(), rel_err_1, time1);

  // TRILINOS PART:

  // Project the initial condition to obtain the initial
  // coefficient vector.
  info("Projecting to obtain initial vector for the Newton's method.");
  CustomInitialSolution sln_tmp(&mesh);
  OGProjection::project_global(&space, &sln_tmp, coeff_vec, matrix_solver);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  CustomWeakForm wf2(JFNK, PRECOND == 1, PRECOND == 2);

  // Initialize DiscreteProblem.
  DiscreteProblem dp2(&wf2, &space);

  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver nox_solver(&dp2, message_type, "GMRES", "Newton", ls_tolerance, "", flag_absresid, abs_resid, 
                       flag_relresid, rel_resid, max_iters);
  nox_solver.set_init_sln(coeff_vec);

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the nonlinear problem using NOX.
  info("Assembling by DiscreteProblem, solving by NOX.");
  Solution sln2;
  if (nox_solver.solve())
  {
    Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln2);
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
         nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
         nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else
    error("NOX failed.");

  // CPU time needed by NOX.
  double time2 = cpu_time.tick().last();

  // Calculate error.
  double rel_err_2 = hermes2d.calc_rel_error(&sln2, &ex, HERMES_H1_NORM) * 100;
  info("Solution 2 (NOX): exact H1 error: %g%% (time %g + %g = %g [s])", rel_err_2, proj_time, time2, proj_time+time2);

  // Show NOX solution.
  ScalarView view2("Solution 2", new WinGeom(510, 0, 500, 400));
  view2.show(&sln2);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
