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

const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.
const int P_INIT = 4;                             // Initial polynomial degree of all mesh elements.

const bool JFNK = true;                          // true = jacobian-free method,
// false = Newton.
const int PRECOND = 1;                            // Preconditioning by jacobian (1) or approximation of jacobian (2)
// in case of JFNK,
// Default ML proconditioner in case of Newton.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
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

double ls_tolerance = 1e-2;                       // Tolerance for linear system.
unsigned flag_absresid = 1;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-8;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 0;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-2;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  MeshFunctionSharedPtr<double> sln(new Hermes::Hermes2D::Solution<double>());

  MeshFunctionSharedPtr<double> ex(new CustomExactSolution(mesh));

  // Initialize the weak formulation for Trilinos.
  WeakFormSharedPtr<double> wf(new CustomWeakForm(JFNK, PRECOND == 1, PRECOND == 2));

  // Initialize the NOX solver with the vector "coeff_vec".
  NewtonSolverNOX<double> nox_solver(wf, space);
  nox_solver.set_verbose_output(true);
  nox_solver.set_output_flags(message_type);
  nox_solver.set_ls_tolerance(ls_tolerance);
  nox_solver.set_conv_abs_resid(abs_resid);
  nox_solver.set_conv_iters(max_iters);

  // Choose preconditioning.
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the nonlinear problem using NOX.
  try
  {
    nox_solver.solve(ex);
  }
  catch (Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }
  Solution<double>::vector_to_solution(nox_solver.get_sln_vector(), space, sln);

  // Show NOX solution.
  Views::ScalarView view("Solution 2", new Views::WinGeom(510, 0, 500, 400));
  view.show(sln);

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}