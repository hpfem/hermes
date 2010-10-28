#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This example solves a system of two linear second-order equations 
// u' + k^2 v = 0
// u - v' = 0
// which is equivalent to u'' + k^2 u = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions 
// u(0) = 0, v(0) = k
// The exact solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

// General input.
static int NEQ = 2;
int NELEM = 20;               // Number of elements.
double A = 0, B = 2*M_PI;     // Domain end points.
int P_init = 2;               // Initial polynomal degree.
double k = 1.0;               // The constant in the equation.

// Newton's method.
double NEWTON_TOL = 1e-5;
int NEWTON_MAX_ITER = 150;

// Boundary conditions.
double Val_dir_left_0 = 0;
double Val_dir_left_1 = k;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"

int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Mesh *mesh = new Mesh(A, B, NELEM, P_init, NEQ);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_0);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_1);
  info("N_dof = %d", mesh->assign_dofs());

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
  int ndof = mesh->get_num_dofs();

  // Fill vector y using dof and coeffs arrays in elements.
  double *y = new double[ndof];
  solution_to_vector(mesh, y);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1)
  {
    // Assemble the Jacobian matrix and residual vector.
    dp->assemble_matrix_and_vector(mesh, matrix, rhs);

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

    // Calculate the coefficient vector.
    bool solved = solver->solve();
    if (solved) 
    {
      double* solution_vector = new double[ndof];
      solution_vector = solver->get_solution();
      for(int i=0; i<ndof; i++) y[i] += solution_vector[i];
      // No need to deallocate the solution_vector here, it is done later by the call to ~Solver.
      solution_vector = NULL;
    }
    it++;

    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    // Copy coefficients from vector y to elements.
    vector_to_solution(y, mesh);
  }
  
  delete matrix;
  delete rhs;
  delete solver;

  // Plot the solution.
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  info("Done.");
  return 1;
}
