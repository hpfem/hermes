#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example explains how to use Newton boundary conditions. Again,
// a Filter is used to visualize the solution gradient.
//
// PDE: Laplace equation -Laplace u = 0 (this equation describes, among
// many other things, also stationary heat transfer in a homogeneous linear
// material).
//
// BC: u = T1 ... fixed temperature on Gamma_left (Dirichlet)
//     du/dn = 0 ... insulated wall on Gamma_outer and Gamma_inner (Neumann)
//     du/dn = H*(u - T0) ... heat flux on Gamma_bottom (Newton)
//
// Note that the last BC can be written in the form  du/dn - H*u = -H*T0.
//
// The following parameters can be changed:

const int P_INIT = 6;                             // Uniform polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int CORNER_REF_LEVEL = 12;                  // Number of mesh refinements towards the re-entrant corner.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_OUTER = 2, BDY_LEFT = 3, BDY_INNER = 4;

// Problem parameters.
const double T1 = 30.0;       // Prescribed temperature on Gamma_left.
const double T0 = 20.0;       // Outer temperature on Gamma_bottom.
const double H  = 0.05;       // Heat flux on Gamma_bottom.

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
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_LEFT);
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_OUTER, BDY_INNER));
  bc_types.add_bc_newton(BDY_BOTTOM);

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_const(BDY_LEFT, T1);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_surf), BDY_BOTTOM);
  wf.add_vector_form_surf(callback(linear_form_surf), BDY_BOTTOM);

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

  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else
    error ("Matrix solver failed.\n");

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&sln);

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", new WinGeom(450, 0, 400, 350));
  MagFilter grad(Hermes::vector<MeshFunction *>(&sln, &sln), 
                 Hermes::vector<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for all views to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}
