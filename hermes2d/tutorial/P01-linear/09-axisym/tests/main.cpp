#include "hermes2d.h"

// This test makes sure that example 06-bc-newton works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

int UNIFORM_REF_LEVEL = 1;                        // number of initial uniform mesh refinements.
int CORNER_REF_LEVEL = 3;                         // number of mesh refinements towards the re-entrant corner.
int P_INIT = 6;                                   // Uniform polynomial degree of all mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double T1 = 300.0;      // Prescribed temperature on Gamma_bottom.
const double T0 = 20.0;       // Outer temperature.
const double LAMBDA = 100;    // Thermal conductivity.
const double h  = 10.0;       // Heat flux on Gamma_heat_flux.

// Boundary markers.
const std::string BDY_HEAT_FLUX = "Heat flux", BDY_BOTTOM = "Bottom", BDY_SYM = "Symmetry";

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY_BOTTOM, T1);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf(h, T0, LAMBDA, BDY_HEAT_FLUX);

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
    if (p_init == 1 && fabs(sum - 1146.15) > 1e-1) success = 0;
    if (p_init == 2 && fabs(sum - 1145.97) > 1e-1) success = 0;
    if (p_init == 3 && fabs(sum - 1145.97) > 1e-1) success = 0;
    if (p_init == 4 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 5 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 6 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 7 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 8 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 9 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 10 && fabs(sum - 1145.96) > 1e-1) success = 0;
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
