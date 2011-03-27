#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This example solves a nonlinear system of four first-order equations.
//
//  PDE: x1' - DAMPING*(1 - x2^2)*x1   +   x2 = 0
//       x2'         -x1   +   x3 = 0
//       x3'         -x2   +   x4 = 0
//       x4'         -x3          = 0.
//
//  Interval: (0, 10).
//
//  BC: Dirichlet, x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0.
//
//  Exact solution: u(x) = exp(x), v(x) = exp(-x).
//
//  The following parameters can be changed:
const int NEQ = 4;                       // Number of equations.
const int NELEM = 500;                   // Number of elements.
const double A = 0, B = 10;              // Domain end points.
const int P_INIT = 2;                    // Polynomial degree.

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
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
Hermes::vector<BCSpec *>DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0,1), new BCSpec(1,0), new BCSpec(2,0), new BCSpec(3,0));
Hermes::vector<BCSpec *> DIR_BC_RIGHT = Hermes::vector<BCSpec *>();

// Weak forms for Jacobi matrix and residual.
#include "definitions.cpp"


int main() 
{
  // Create space, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ);

  // Enumerate basis functions, info for user.
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf(4);
  wf.add_matrix_form(0, 0, jacobian_1_1);
  wf.add_matrix_form(0, 1, jacobian_1_2);
  wf.add_matrix_form(1, 0, jacobian_2_1);
  wf.add_matrix_form(1, 1, jacobian_2_2);
  wf.add_matrix_form(1, 2, jacobian_2_3);
  wf.add_matrix_form(2, 1, jacobian_3_2);
  wf.add_matrix_form(2, 2, jacobian_3_3);
  wf.add_matrix_form(2, 3, jacobian_3_4);
  wf.add_matrix_form(3, 2, jacobian_4_3);
  wf.add_matrix_form(3, 3, jacobian_4_4);
  wf.add_vector_form(0, residual_1);
  wf.add_vector_form(1, residual_2);
  wf.add_vector_form(2, residual_3);
  wf.add_vector_form(3, residual_4);

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

  // Cleanup.
  for(unsigned j = 0; j < DIR_BC_LEFT.size(); j ++)
      delete DIR_BC_LEFT[j];
  DIR_BC_LEFT.clear();

  delete matrix;
  delete rhs;
  delete solver;
  delete[] coeff_vec;
  delete space;

  info("Done.");
  return 0;
}
