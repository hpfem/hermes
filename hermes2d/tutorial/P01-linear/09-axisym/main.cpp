#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to solve exisymmetric problems. The domain of interest
// is a hollow cylinder whose axis is aligned with the y-axis. It has fixed
// temperature on the bottom face, and a Newton-type heat flux condition on
// all other faces.
//
// PDE: Laplace equation -Laplace u = 0 (y-axis is the axis of symmetry).
//
// BC: u = T1 ... fixed temperature on Gamma_bottom,
//     -LAMBDA * du/dn = ALPHA*(u - T0) ... heat flux on the rest of the boundary (Gamma_heat_flux).
//
// Note that the last BC can be written in the form  du/dn - H*u = -H*T0.
//
// The following parameters can be changed:

const int P_INIT = 4;                             // Uniform polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Problem parameters.
const double T1 = 30.0;       // Prescribed temperature on Gamma_bottom.
const double T0 = 20.0;       // Outer temperature.
const double LAMBDA = 386;    // Thermal conductivity.
const double ALPHA = 5.0;     // Heat flux coefficient on Gamma_heat_flux.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential("Bottom", T1);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf(LAMBDA, ALPHA, T0, "Heat flux");

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution sln;
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 300, 400));
  view.show(&sln);
  ScalarView gradview("Gradient", new WinGeom(310, 0, 300, 400));
  MagFilter grad(Hermes::vector<MeshFunction *>(&sln, &sln), 
                 Hermes::vector<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for all views to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  delete [] coeff_vec;

  return 0;
}
