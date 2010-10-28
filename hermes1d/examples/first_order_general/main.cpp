#define HERMES_REPORT_ALL
#include "hermes1d.h"

//  This example solves a general first-order equation 
//  y' = f(y, x) in an interval (A, B), equipped with the 
//  initial condition y(A) = YA. The function f can be linear
//  or nonlinear in 'y', as long as it is differentiable
//  with respect to this variable (needed for the Newton's method). 
//
//  Equation: y' = f(y, x)
//
//  Initial condition: y(A) = YA
//
//  Interval: (A, B)
//
//  The following parameters can be changed.

int P_INIT = 2;                                   // Initial polynomal degree.
static int NEQ = 1;                               // Number of equations.
int NELEM = 10;                                   // Number of elements.
double A = 0, B = 10;                             // Domain end points.
double YA = 1;                                    // Equation parameter.
double NEWTON_TOL = 1e-5;                         // Tolerance for the Newton's method.
int NEWTON_MAX_ITER = 150;                        // Maximum allowed number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Function f(y, x).
double f(double y, double x) {
  return -y;
}

// Function dfdy(y, x).
double dfdy(double y, double x) {
  return -1;
}

// Weak forms for the Jacobi matrix and residual.
#include "forms.cpp"

int main() {
  // Create mesh, set Dirichlet BC, enumerate basis functions.
  Mesh *mesh = new Mesh(A, B, NELEM, P_INIT, NEQ);
  mesh->set_bc_left_dirichlet(0, YA);
  info("N_dof = %d.", mesh->assign_dofs());

  // Initialize the weak formulation.
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Obtain the number of degrees of freedom.
  int ndof = mesh->get_num_dofs();

  // Copy solution coefficients into a vector.
  double *y = new double[ndof];
  solution_to_vector(mesh, y);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Newton's iteration.
  int it = 1;
  while (1)
  {
    // Assemble the stiffness matrix and right-hand side vector.
    dp->assemble_matrix_and_vector(mesh, matrix, rhs);

    // Calculate L2 norm of the residual vector.
    double res_norm_squared = 0;
    for(int i = 0; i < ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);
    double res_norm = sqrt(res_norm_squared);

    info("---- Newton iter %d, residual norm: %.15f.", it, res_norm);

    // If l2 norm of the residual vector is within tolerance, then quit.
    // Latest solution is in the vector 'y'.
    if(res_norm < NEWTON_TOL && it > 1) break;

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
  
  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  // Plot the solution.
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  // Plot the resulting mesh.
  mesh->plot("mesh.gp");

  info("Done.");
  return 1;
}
