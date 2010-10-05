#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

using namespace RefinementSelectors;

//  This test makes sure that example 16-newton-elliptic-2 works correctly.

const int P_INIT = 2;                             // Initial polynomial degree
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 8;                  // Maximum allowed number of Newton iterations
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  // Using the Dirichlet lift elevated by two
  double val = dir_lift(x, y, dx, dy) + 2;
  return val;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation
  WeakForm wf;
  wf.add_matrix_form(callback(jac), HERMES_UNSYM, HERMES_ANY);
  wf.add_vector_form(callback(res), HERMES_ANY);

  // Initialize the FE problem.
  bool is_linear = false;
  FeProblem fep(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[get_num_dofs(&space)] ;
  Solution* init_sln = new Solution(&mesh, init_cond);
  project_global(&space, init_sln, coeff_vec, matrix_solver); 
  delete init_sln;

  // Perform Newton's iteration.
  int it = 1;
  bool success = false;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = get_num_dofs(&space);

    // Assemble the Jacobian matrix and residual vector.
    fep.assemble(coeff_vec, matrix, rhs, false);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
    
    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, get_num_dofs(&space), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

    // Solve the linear system.
    if(!(success = solver->solve()))
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
    
    if (it >= NEWTON_MAX_ITER)
      error ("Newton method did not converge.");

    it++;
  }

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Cleanup.
  delete []coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success) {  // should pass with NEWTON_MAX_ITER = 8 and fail with NEWTON_MAX_ITER = 7
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

