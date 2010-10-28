#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This example solves the mathematical pendulum equation 
// y'' + k**2 * sin(y) = 0 in an interval (A, B), equipped with the 
// initial conditions y(A) = Init_angle, y'(0) = Init_vel. The 
// system is decomposed into two first order ODE and solved via 
// the Newton's method starting from zero initial condition.
// Note that the method diverges for longer time intervals, 
// depending on the interval length, number of elements, and 
// the initial polynomial degree.
//
// Derivation:
// m*l*u'' = -m*g*sin(u)
// so:
// u'' + k^2 * sin(u) = 0
// with k^2 = g/l
// so we have to solve a system of two nonlinear second-order equations
// v' + k^2 sin u = 0
// u' - v = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions
// u(0) = 0, v(0) = k
// The approximate (linearized) solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

// General input.
static int NEQ = 2;
int NELEM = 1292;             // Number of elements.
double A = 0, B = 10;         // Domain end points.
int P_init = 1;               // Initial polynomal degree.
double k = 0.5;

// Newton's method.
double NEWTON_TOL = 1e-5;
int NEWTON_MAX_ITER = 150;

// Boundary conditions.
double Init_angle = M_PI/2.;  // Initial angle.
double Init_vel = 0;          // Initial velocity.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"


int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space *space = new Space(A, B, NELEM, P_init, NEQ);
  space->set_bc_left_dirichlet(0, Init_angle);
  space->set_bc_left_dirichlet(1, Init_vel);
  info("N_dof = %d", space->assign_dofs());

  // Initialize the FE problem.
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 1, jacobian_0_1);
  dp->add_matrix_form(1, 0, jacobian_1_0);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);

  // Newton's loop.
  // Obtain the number of degrees of freedom.
  int ndof = Space::get_num_dofs(space);

  // Fill vector y using dof and coeffs arrays in elements.
  double *coeff_vec = new double[ndof];
  solution_to_vector(space, coeff_vec);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1)
  {
    // Assemble the Jacobian matrix and residual vector.
    dp->assemble_matrix_and_vector(space, matrix, rhs);

    // Calculate the l2-norm of residual vector.
    double res_norm_squared = 0;
    for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

    // Info for user.
    info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_norm_squared < NEWTON_TOL*NEWTON_TOL && it > 1) break;

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
    
    // Copy coefficients from vector y to elements.
    vector_to_solution(coeff_vec, space);

    it++;
  }
  
  // Cleanup.
  delete matrix;
  delete rhs;
  delete solver;

  // Plot the solution.
  Linearizer l(space);
  l.plot_solution("solution.gp");

  info("Done.");
  return 1;
}
