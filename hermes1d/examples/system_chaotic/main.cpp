#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This example solves a nonlinear system of four first-order equations
// x1' - DAMPING*(1 - x2^2)*x1   +   x2 = 0
// x2'         -x1   +   x3 = 0
// x3'         -x2   +   x4 = 0
// x4'         -x3          = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

// General input.
static int NEQ = 4;
int NELEM = 500;            // Number of elements.
double A = 0, B = 10;       // Domain end points.
int P_init = 2;             // Initial polynomal degree.

// Damping parameter.
int DAMPING_STEPS = 20;     // Number of damping steps. The entire problem
                            // will be run repeatedly, with the DAMPING parameter 
                            // increased from 0 to 1 in DAMPING_STEPS. Every time, 
                            // the last result is used as initial cond. for the 
                            // new computation.   
double DAMPING = 1.0;       // DAMPING is an artificial param. used to 
                            // reduce the strength of the nonlinearity. 
                            // (The nonlinearity is multiplied with it.)

// Newton's method.
double NEWTON_TOL = 1e-5;
int NEWTON_MAX_ITER = 150;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Boundary conditions.
double Val_dir_left_1 = 1;
double Val_dir_left_2 = 0;
double Val_dir_left_3 = 0;
double Val_dir_left_4 = 0;

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"


int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space *space = new Space(A, B, NELEM, P_init, NEQ);
  space->set_bc_left_dirichlet(0, Val_dir_left_1);
  space->set_bc_left_dirichlet(1, Val_dir_left_2);
  space->set_bc_left_dirichlet(2, Val_dir_left_3);
  space->set_bc_left_dirichlet(3, Val_dir_left_4);
  info("N_dof = %d", space->assign_dofs());

  // Initialize the FE problem.
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_1_1);
  dp->add_matrix_form(0, 1, jacobian_1_2);
  dp->add_matrix_form(1, 0, jacobian_2_1);
  dp->add_matrix_form(1, 1, jacobian_2_2);
  dp->add_matrix_form(1, 2, jacobian_2_3);
  dp->add_matrix_form(2, 1, jacobian_3_2);
  dp->add_matrix_form(2, 2, jacobian_3_3);
  dp->add_matrix_form(2, 3, jacobian_3_4);
  dp->add_matrix_form(3, 2, jacobian_4_3);
  dp->add_matrix_form(3, 3, jacobian_4_4);
  dp->add_vector_form(0, residual_1);
  dp->add_vector_form(1, residual_2);
  dp->add_vector_form(2, residual_3);
  dp->add_vector_form(3, residual_4);

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
