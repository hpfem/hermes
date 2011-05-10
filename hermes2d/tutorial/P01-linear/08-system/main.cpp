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
// BC: du_1/dn = f0 on Gamma_top (top edge),
//     du_2/dn = f1 on Gamma_top (top edge),
//     u_1 = 0 and u_2 = 0 on Gamma_bottom (bottom edge),
//     du_1/dn = 0 on Gamma_rest (rest of boundary),
//     du_2/dn = 0 on Gamma_rest (rest of boundary).
//
// The following parameters can be changed:

const int P_INIT = 6;                                      // Initial polynomial degree of all elements.
const double NEWTON_TOL = 1e-8;                            // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                           // Maximum allowed number of Newton iterations.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                           // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double rho = 8000.0;                                 // Density.
const double g1 = -9.81;                                   // Gravitational acceleration.
const double f0  = 0;                                      // Surface force in x-direction.
const double f1  = 8e4;                                    // Surface force in y-direction.

// Weak forms.
#include "definitions.h"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh, mesh1;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform uniform mesh refinement.
  mesh.refine_all_elements();

  // Show mesh.
  MeshView mv("Mesh", new WinGeom(0, 0, 580, 400));
  mv.show(&mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_disp("Bottom", 0.0);
  EssentialBCs bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u1_space(&mesh, &bcs, P_INIT);
  H1Space u2_space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space));
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "Top", f0, f1);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u1_space, &u2_space));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  bool verbose = true;
  bool jacobian_changed = true;
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution u1_sln, u2_sln;
  Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&u1_space, &u2_space), 
                                Hermes::vector<Solution *>(&u1_sln, &u2_sln));
  
  // Visualize the solution.
  ScalarView view("Von Mises stress [Pa]", new WinGeom(590, 0, 700, 400));
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // First Lame constant.
  double mu = E / (2*(1 + nu));                        // Second Lame constant.
  VonMisesFilter stress(Hermes::vector<MeshFunction *>(&u1_sln, &u2_sln), lambda, mu);
  view.show_mesh(false);
  view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.5e5);

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete [] coeff_vec;
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

