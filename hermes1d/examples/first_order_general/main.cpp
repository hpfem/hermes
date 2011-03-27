#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

//  This example solves a general first-order equation 
//  y' = f(y, x) in an interval (A, B), equipped with the 
//  initial condition y(A) = YA. The function f can be linear
//  or nonlinear in 'y', as long as it is differentiable
//  with respect to this variable (needed for the Newton's method).
//
//  PDE: y' = f(y, x).
//  The function f(y, x) can be changed below.
//
//  Interval: (A, B).
//
//  IC: y(A) = YA.
//
//  DC: Determined by the initial condition.
//
//  The following parameters can be changed:

const int NEQ = 1;                                // Number of equations.
const int NELEM = 10;                             // Number of elements.
const double A = 0, B = 10;                       // Domain end points.
const double YA = 1;                              // Equation parameter.
const int P_INIT = 2;                             // Polynomial degree.
double NEWTON_TOL = 1e-5;                         // Tolerance.
int NEWTON_MAX_ITER = 150;                        // Max. number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
BCSpec DIR_BC_LEFT(0, YA);
BCSpec DIR_BC_RIGHT;

// Function f(y, x).
double f(double y, double x) {
  return -y*y;
}

// Function dfdy(y, x).
double dfdy(double y, double x) {
  return -2*y;
}

// Weak forms for the Jacobi matrix and residual.
#include "definitions.cpp"

int main() 
{
  // Create space, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ);
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jacobian);
  wf.add_vector_form(residual);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);

  // Set zero initial condition.
  double *coeff_vec = new double[ndof];
  set_zero(coeff_vec, ndof);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1) 
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(space);

    // Assemble the Jacobian matrix and residual vector.
    dp->assemble(coeff_vec, matrix, rhs);

    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(space), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_l2_norm < NEWTON_TOL && it > 1) break;

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

    // Solve the linear system.
    if(!solver->solve())
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];

    // If the maximum number of iteration has been reached, then quit.
    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    it++;
  }

  // Plot the solution.
  Linearizer l(space);
  l.plot_solution("solution.dat");

  // cleaning
  delete dp;
  delete rhs;
  delete solver;
  delete[] coeff_vec;
  delete space;
  delete matrix;

  info("Done.");
  return 0;
}
