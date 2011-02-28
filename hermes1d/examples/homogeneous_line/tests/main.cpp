#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that example "homogeneous_line" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

double L=25e-9;               // induktance [H/m]
double C=10e-12;              // capacitance [F/m]
double G=1e-9;                // conductance [S/m]
double R=1e-3;                // resistance [Ohm/m]
double l=10 ;                 // length of the line [m]
double omega=2*M_PI*3e8;      //
double Zl=60;                 // load impedance[Ohm]

//double Val_newton_alpha_U_Re=-R/Zl;
//double Val_newton_alpha_U_Im=-omega*L/Zl;
//double Val_newton_alpha_I_Re=-G*Zl;
//double Val_newton_alpha_I_Im=-omega*C*Zl;

// General input:
static int NEQ = 4;
int NELEM = 1000;          // number of elements
double A = 0, B = l;        // domain end points
int P_INIT = 2;             // initial polynomial degree

// Matrix solver.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  // Only relevant for iterative matrix solvers:

// Newton's method.
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.

// Boundary conditions.
Hermes::vector<BCSpec *> DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0,1), new BCSpec(0,0), 
                              new BCSpec(0,0), new BCSpec(0,0));
Hermes::vector<BCSpec *> DIR_BC_RIGHT = Hermes::vector<BCSpec *>();

//At the end of the line is an indirect boundary condition U(l) = I(l)*Zl see below

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

int main() 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Create space, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ);
  info("ndof = %d", Space::get_num_dofs(space));

  // Initialize the weak formulation.
  WeakForm wf(4);
  wf.add_matrix_form(0, 0, jacobian_1_1);
  wf.add_matrix_form(0, 2, jacobian_1_3);
  wf.add_matrix_form(0, 3, jacobian_1_4);
  wf.add_matrix_form(1, 1, jacobian_2_2);
  wf.add_matrix_form(1, 2, jacobian_2_3);
  wf.add_matrix_form(1, 3, jacobian_2_4);
  wf.add_matrix_form(2, 0, jacobian_3_1);
  wf.add_matrix_form(2, 1, jacobian_3_2);
  wf.add_matrix_form(2, 2, jacobian_3_3);
  wf.add_matrix_form(3, 0, jacobian_4_1);
  wf.add_matrix_form(3, 1, jacobian_4_2);
  wf.add_matrix_form(3, 3, jacobian_4_4);
  wf.add_vector_form(0, residual_1);
  wf.add_vector_form(1, residual_2);
  wf.add_vector_form(2, residual_3);
  wf.add_vector_form(3, residual_4);
  wf.add_matrix_form_surf(0, 0, jacobian_surf_right_U_Re, BOUNDARY_RIGHT);
  wf.add_matrix_form_surf(0, 2, jacobian_surf_right_U_Im, BOUNDARY_RIGHT);
  wf.add_matrix_form_surf(1, 1, jacobian_surf_right_I_Re, BOUNDARY_RIGHT);
  wf.add_matrix_form_surf(1, 3, jacobian_surf_right_I_Im, BOUNDARY_RIGHT);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);

  // Set zero initial condition.
  double *coeff_vec = new double[Space::get_num_dofs(space)];
  set_zero(coeff_vec, Space::get_num_dofs(space));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  bool success = false;
  while (1) {
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
    for(int i = 0; i < Space::get_num_dofs(space); i++) rhs->set(i, -rhs->get(i));

    // Solve the linear system.
    if(!(success = solver->solve()))
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < Space::get_num_dofs(space); i++) coeff_vec[i] += solver->get_solution()[i];

    // If the maximum number of iteration has been reached, then quit.
    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    it++;
  }

  info("Total running time: %g s", cpu_time.accumulated());

  // Test variable.
  info("ndof = %d.", Space::get_num_dofs(space));

  if (success)
  {
    info("Success!");
    return ERROR_SUCCESS;
  }
  else
  {
    info("Failure!");
    return ERROR_FAILURE;
  }
}
