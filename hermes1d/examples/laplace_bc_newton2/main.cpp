#define HERMES_REPORT_ALL
#include "hermes1d.h"

//  This example solves the Poisson equation.
//
//  PDE: -u'' - f = 0.
//
//  Interval: (A, B).
//
//  BC: Newton.
//
//  The following parameters can be changed:
const int NEQ = 1;                      // Number of equations.
const int NELEM = 3;                    // Number of elements.
const double A = 0, B = 2*M_PI;         // Domain end points.
const int P_INIT = 3;                   // Polynomial degree.

// Newton's method.
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.


MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
double Val_newton_alpha_left = 2;
double Val_newton_beta_left = -2;
double Val_newton_alpha_right = 1;
double Val_newton_beta_right = 1;

// Function f(x).
double f(double x) {
  return sin(x);
}

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"

int main() 
{
  // Create space, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, P_INIT, NEQ);
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jacobian_vol);
  wf.add_vector_form(residual_vol);
  wf.add_matrix_form_surf(jacobian_surf_left, BOUNDARY_LEFT);
  wf.add_vector_form_surf(residual_surf_left, BOUNDARY_LEFT);
  wf.add_matrix_form_surf(jacobian_surf_right, BOUNDARY_RIGHT);
  wf.add_vector_form_surf(residual_surf_right, BOUNDARY_RIGHT);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);

  // Newton's loop.
  // Fill vector coeff_vec using dof and coeffs arrays in elements.
  double *coeff_vec = new double[Space::get_num_dofs(space)];
  get_coeff_vector(space, coeff_vec);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1) {
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
    
    // Copy coefficients from vector y to elements.
    set_coeff_vector(coeff_vec, space);

    it++;
  }
  
  // Plot the solution.
  Linearizer l(space);
  l.plot_solution("solution.gp");

  // Plot the resulting space.
  space->plot("space.gp");

  info("Done.");
  return 0;
}
