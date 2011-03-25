#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example explains how to create multiple spaces over a mesh and use them
// to solve a simple problem of linear elasticity. Note how Tuples are used, 
// they replace variable-length argument lists. At the end, VonMises filter is 
// used to visualize the stress.
//
// PDE: Lame equations of linear elasticity.
//
// BC: du_1/dn = f0 on Gamma_3 (top edge),
//     du_2/dn = f1 on Gamma_3 (top edge),
//     u_1 = 0 and u_2 = 0 on Gamma_1 (bottom edge),
//     du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5 (rest of boundary),
//     du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5 (rest of boundary).
//
// The following parameters can be changed:

const int P_INIT = 6;                                      // Initial polynomial degree of all elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                           // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_1 = "1", BDY_3 = "3";

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double rho = 8000.0;                                 // Density.
const double g1 = -9.81;                                   // Gravitational acceleration.
const double f0  = 0;                                      // Surface force in x-direction.
const double f1  = 8e4;                                    // Surface force in y-direction.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, mesh1;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform uniform mesh refinement.
  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_disp(BDY_1, 0.0);
  EssentialBCs bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u1_space(&mesh, &bcs, P_INIT);
  H1Space u2_space(&mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)));

  // Initialize the weak formulation.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, BDY_3, f0, f1);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u1_space, &u2_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solutions.
  Solution u1_sln, u2_sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solutions.
  info("Solving the matrix problem.");
  if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(&u1_space, &u2_space), 
                                                    Hermes::vector<Solution *>(&u1_sln, &u2_sln));
  else error ("Matrix solver failed.\n");
  
  // Visualize the solution.
  ScalarView view("Von Mises stress [Pa]", new WinGeom(0, 0, 800, 400));
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // First Lame constant.
  double mu = E / (2*(1 + nu));                        // Second Lame constant.
  VonMisesFilter stress(Hermes::vector<MeshFunction *>(&u1_sln, &u2_sln), lambda, mu);
  view.show_mesh(false);
  view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.5e5);

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

