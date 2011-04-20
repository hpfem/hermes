#include "hermes2d.h"

// This test makes sure that example 06-bc-newton works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

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
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential("Bottom", T1);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf(LAMBDA, ALPHA, T0, "Heat flux");

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  for (int p_init = 1; p_init <= 10; p_init++) {

    printf("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);
    int ndof = Space::get_num_dofs(&space);
    info("ndof = %d", ndof);

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

    double sum = 0;
    for (int i=0; i < ndof; i++) sum += coeff_vec[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 597.764) > 1e-1) success = 0;
    if (p_init == 2 && fabs(sum - 597.808) > 1e-1) success = 0;
    if (p_init == 3 && fabs(sum - 597.8) > 1e-1) success = 0;
    if (p_init == 4 && fabs(sum - 597.807) > 1e-1) success = 0;
    if (p_init == 5 && fabs(sum - 597.8) > 1e-1) success = 0;
    if (p_init == 6 && fabs(sum - 597.806) > 1e-1) success = 0;
    if (p_init == 7 && fabs(sum - 597.8) > 1e-1) success = 0;
    if (p_init == 8 && fabs(sum - 597.806) > 1e-1) success = 0;
    if (p_init == 9 && fabs(sum - 597.8) > 1e-1) success = 0;
    if (p_init == 10 && fabs(sum - 597.806) > 1e-1) success = 0;

    delete [] coeff_vec;
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
