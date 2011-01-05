#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

using namespace RefinementSelectors;

//  This test makes sure that example 15-picard works correctly.

const int P_INIT = 2;                             // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;                   // Number of initial refinements towards boundary.
const double PICARD_TOL = 1e-6;                   // Stopping criterion for the Picard's method.
const int MAX_PICARD_ITER_NUM = 100;              // Maximum allowed number of Picard iterations.
const double INIT_COND_CONST = 3.0;               // Constant initial condition.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u) 
{ 
  return 1 + pow(u, 4); 
}

// Right-hand side (can be a general function of 'x' and 'y').
template<typename Real>
Real rhs(Real x, Real y)
{
  return 1.0;
}

// Boundary markers.
const int BDY_DIRICHLET = 1;

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

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Declare solutions.
  Solution sln, *sln_prev;

  // Initialize previous iteration solution for the Picard method.
  sln_prev = new Solution();
  sln_prev->set_const(&mesh, INIT_COND_CONST);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), HERMES_NONSYM, HERMES_ANY, sln_prev);
  wf.add_vector_form(callback(linear_form), HERMES_ANY);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Picard's iteration.
  double rel_err;
  int iter_count = 0;
  while (true) {
    // Assemble the stiffness matrix and right-hand side.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else error ("Matrix solver failed.\n");

    double rel_error = calc_rel_error(sln_prev, &sln, HERMES_H1_NORM) * 100;
    info("Relative error: %g%%", rel_error);

    // Stopping criterion.
    if (rel_error < PICARD_TOL || iter_count >= MAX_PICARD_ITER_NUM) break;

    // Saving solution for the next iteration;
    sln_prev->copy(&sln);
   
    iter_count++;
  }

  if (iter_count >= MAX_PICARD_ITER_NUM) warn("Picard's iteration did not converge.");
  else info("Picard's iteration needed %d steps.", iter_count);

  // Cleanup.
  delete matrix;
  delete rhs;
  delete solver;

  ndof = Space::get_num_dofs(&space);
  info("Coordinate (-10, -10) value = %lf", sln.get_pt_value(-10.0, -10.0));
  info("Coordinate ( -6,  -6) value = %lf", sln.get_pt_value(-6.0, -6.0));
  info("Coordinate ( -2,  -2) value = %lf", sln.get_pt_value(-2.0, -2.0));
  info("Coordinate (  2,   2) value = %lf", sln.get_pt_value(2.0, 2.0));
  info("Coordinate (  6,   6) value = %lf", sln.get_pt_value(6.0, 6.0));
  info("Coordinate ( 10,  10) value = %lf", sln.get_pt_value(10.0, 10.0));

  double coor_x_y[6] = {-10.0, -6.0, -2.0, 2.0, 6.0, 10.0};
  double value[6] = {0.000000, 2.255137, 2.623581, 2.623581, 2.255137, 0.000000};
  for (int i = 0; i < 6; i++)
  {
    if ((value[i] - sln.get_pt_value(coor_x_y[i], coor_x_y[i])) < 1E-6)
    {
      printf("Success!\n");
    }
    else
    {
      printf("Failure!\n");
      return ERR_FAILURE;
    }
  }
  return ERR_SUCCESS;
}

