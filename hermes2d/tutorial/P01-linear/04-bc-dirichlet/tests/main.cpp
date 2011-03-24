#include "hermes2d.h"

// This test makes sure that example 04-bc-dirichlet works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

int P_INIT = 2;                                   // Initial polynomial degree in all elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double CONST_F = -4.0;                            // constant right-hand side.

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);
  mesh.refine_all_elements();

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set exact solution.
  CustomExactSolution exact(&mesh, CONST_F);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf(CONST_F);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst bc_essential("Dirichlet", &exact);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  Solution sln;
  for (int p_init = 1; p_init <= 10; p_init++) {

    printf("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);

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
    bool rhsonly = false;
    dp.assemble(matrix, rhs, rhsonly);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else
      error ("Matrix solver failed.\n");

    int ndof = Space::get_num_dofs(&space);
    printf("ndof = %d\n", ndof);
    double sum = 0;
    for (int i=0; i < ndof; i++) sum += solver->get_solution()[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 1.7251) > 1e-3) success = 0;
    if (p_init == 2 && fabs(sum - 3.79195) > 1e-3) success = 0;
    if (p_init == 3 && fabs(sum - 3.80206) > 1e-3) success = 0;
    if (p_init == 4 && fabs(sum - 3.80156) > 1e-3) success = 0;
    if (p_init == 5 && fabs(sum - 3.80155) > 1e-3) success = 0;
    if (p_init == 6 && fabs(sum - 3.80154) > 1e-3) success = 0;
    if (p_init == 7 && fabs(sum - 3.80154) > 1e-3) success = 0;
    if (p_init == 8 && fabs(sum - 3.80153) > 1e-3) success = 0;
    if (p_init == 9 && fabs(sum - 3.80152) > 1e-3) success = 0;
    if (p_init == 10 && fabs(sum - 3.80152) > 1e-3) success = 0;
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
