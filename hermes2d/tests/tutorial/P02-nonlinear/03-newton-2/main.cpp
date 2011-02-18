#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This test makes sure that example 17-newton-elliptic-2 works correctly.

const int P_INIT = 2;                             // Initial polynomial degree
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 9;                  // Maximum allowed number of Newton iterations
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

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

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y)
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
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize the weak formulation
  WeakForm wf;
  wf.add_matrix_form(callback(jac), HERMES_NONSYM, HERMES_ANY);
  wf.add_vector_form(callback(res), HERMES_ANY);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)] ;
  Solution* init_sln = new Solution(&mesh, init_cond);
  OGProjection::project_global(&space, init_sln, coeff_vec, matrix_solver); 
  delete init_sln;

  // Perform Newton's iteration.
  bool verbose = true;
  if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Cleanup.
  delete []coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
  info("ndof = %d", ndof);
  info("Coordinate (1, 0) value = %lf", sln.get_pt_value(1.0, 0.0));
  info("Coordinate (3, 0) value = %lf", sln.get_pt_value(3.0, 0.0));
  info("Coordinate (5, 0) value = %lf", sln.get_pt_value(5.0, 0.0));
  info("Coordinate (7, 0) value = %lf", sln.get_pt_value(7.0, 0.0));

  double coor_x[4] = {1.0, 3.0, 5.0, 7.0};
  double coor_y = 0.0;
  double t_value[4] = {2.867436, 2.873677, 2.832594, 2.709390};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (abs(t_value[i] - sln.get_pt_value(coor_x[i], coor_y)) > 1E-6) success = false;
  }

  if (success) {  // should pass with NEWTON_MAX_ITER = 9 and fail with NEWTON_MAX_ITER = 8
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }


}

