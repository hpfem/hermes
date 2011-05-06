#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This test makes sure that example 08-system works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).
const int P_INIT = 8;                                      // Initial polynomial degree of all elements.
const double NEWTON_TOL = 1e-8;                            // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                           // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                           // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double rho = 8000.0;                                 // Density.
const double g1 = -9.81;                                   // Gravitational acceleration.
const double f0  = 0;                                      // External force in x-direction.
const double f1  = 8e4;                                    // External force in y-direction.

// Weak forms.
#include "../definitions.h"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform uniform mesh refinement.
  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_disp("Bottom", 0.0);
  EssentialBCs bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u1_space(&mesh, &bcs, P_INIT);
  H1Space u2_space(&mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)));

  // Initialize the weak formulation.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "Top", f0, f1);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  Solution xsln, ysln;
  for (int p_init = 1; p_init <= 6; p_init++) {
    printf("********* p_init = %d *********\n", p_init);
    u1_space.set_uniform_order(p_init);
    u2_space.set_uniform_order(p_init);
    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space));
    info("ndof = %d", ndof);

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

    double sum = 0;
    for (int i=0; i < ndof; i++) sum += coeff_vec[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 1.41886e-05) > 1e-5) success = 0;
    if (p_init == 2 && fabs(sum - 1.60006e-05) > 1e-5) success = 0;
    if (p_init == 3 && fabs(sum - 1.60810e-05) > 1e-5) success = 0;
    if (p_init == 4 && fabs(sum - 1.61106e-05) > 1e-5) success = 0;
    if (p_init == 5 && fabs(sum - 1.61065e-05) > 1e-5) success = 0;
    if (p_init == 6 && fabs(sum - 1.61112e-05) > 1e-5) success = 0;

    delete [] coeff_vec;
    delete solver;
    delete matrix;
    delete rhs;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

