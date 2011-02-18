#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example illustrates how to use nonhomogeneous (nonzero)
// Dirichlet boundary conditions.
//
// PDE: Poisson equation -Laplace u = CONST_F, where CONST_F is
// a constant right-hand side. It is not difficult to see that
// the function u(x,y) = (-CONST_F/4)*(x^2 + y^2) satisfies the
// above PDE. Since also the Dirichlet boundary conditions
// are chosen to match u(x,y), this function is the exact
// solution.
//
// Note that since the exact solution is a quadratic polynomial,
// Hermes will compute it exactly if all mesh elements are quadratic
// or higher (then the exact solution lies in the finite element space).
// If some elements in the mesh are linear, Hermes will only find
// an approximation, Below you can play with the parameters CONST_F,
// P_INIT, and INIT_REF_NUM.

const int P_INIT = 2;                             // Initial polynomial degree in all elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_OUTER = 2, BDY_LEFT = 3, BDY_INNER = 4;

// Problem parameters.
const double CONST_F = -4.0; 

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y)
{
  return (-CONST_F/4.0)*(x*x + y*y);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_OUTER, BDY_LEFT, BDY_INNER));

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_function(Hermes::vector<int>(BDY_BOTTOM, BDY_OUTER, BDY_LEFT, BDY_INNER), essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the and solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else
    error ("Matrix solver failed.\n");

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&sln);

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}
